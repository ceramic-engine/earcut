package earcut;
 
// Ported from https://github.com/mapbox/earcut by Krtolica Vujadin
// Edited by Jeremy Faivre to optimise allocations and update from latest mapbox earcut fixes

class Earcut {
    
    static public function earcut(data:Array<Float>, ?holeIndices:Array<Int>, dim:Int = 2, ?triangles:Array<Int>):Array<Int> {
        var hasHoles:Bool = holeIndices != null && holeIndices.length > 0;
        var	outerLen:Int = hasHoles ? Std.int(holeIndices[0] * dim) : data.length;
        var	outerNode:Node = linkedList(data, 0, outerLen, dim, true);

        if (triangles == null) triangles = [];

        // Empty triangles data
        if (triangles.length > 0) {
            #if cpp
            untyped triangles.__SetSize(0);
            #else
            triangles.splice(0, triangles.length);
            #end
        }
        
        if (outerNode == null || outerNode.next == outerNode.prev) {
            return triangles;
        }
        
        var minX:Float = 0;
        var minY:Float = 0;
        var maxX:Float = 0;
        var maxY:Float = 0;
        var x:Float = 0;
        var y:Float = 0;
        var invSize:Float = 0;
        
        if (hasHoles) {
            outerNode = eliminateHoles(data, holeIndices, outerNode, dim);
        }
        
        // If the shape is not too simple, we'll use z-order curve hash later; calculate polygon bbox
        if (data.length > 80 * dim) {
            minX = maxX = data[0];
            minY = maxY = data[1];
            
            var i = dim;
            while (i < outerLen) {
                x = data[i];
                y = data[i + 1];
                if (x < minX) {
                    minX = x;
                }
                if (y < minY) {
                    minY = y;
                }
                if (x > maxX) {
                    maxX = x;
                }
                if (y > maxY) {
                    maxY = y;
                }
                
                i += dim;
            }
            
            // minX, minY and invSize are later used to transform coords into integers for z-order calculation
            invSize = Math.max(maxX - minX, maxY - minY);
            invSize = invSize != 0 ? 1 / invSize : 0;
        }
        
        earcutLinked(outerNode, triangles, dim, minX, minY, invSize);

        Node.recycleAllNodes();
        
        return triangles;
    }
    
    /**
     * Create a circular doubly linked list from polygon points in the specified winding order
     */
    static function linkedList(data:Array<Float>, start:Int, end:Int, dim:Int, clockwise:Bool):Node {
        var last:Node = null;
        
        if (clockwise == (signedArea(data, start, end, dim) > 0)) {
            var i:Int = start;
            while (i < end) {
                last = insertNode(i, data[i], data[i + 1], last);
                
                i += dim;
            }
        } 
        else {
            var i:Int = end - dim;
            while (i >= start) {
                last = insertNode(i, data[i], data[i + 1], last);
                i -= dim;
            }
        }
        
        if (last != null && equals(last, last.next)) {
            removeNode(last);
            last = last.next;
        }
        
        return last;
    }

    /**
     * Eliminate colinear or duplicate points
     */
    static function filterPoints(?start:Node, ?end:Node):Node {
        if (start == null) {
            return start;
        }
        if (end == null) {
            end = start;
        }
        
        var p:Node = start;
        var again:Bool = false;
        do {
            again = false;
            
            if (!p.steiner && (equals(p, p.next) || area(p.prev, p, p.next) == 0)) {
                removeNode(p);
                p = end = p.prev;
                if (p == p.next) {
                    break;
                }
                again = true;
            } 
            else {
                p = p.next;
            }
        } 
        while (again || p != end);
        
        return end;
    }

    /**
     * Main ear slicing loop which triangulates a polygon (given as a linked list)
     */
    static function earcutLinked(ear:Node, triangles:Array<Int>, dim:Float, ?minX:Float, ?minY:Float, ?invSize:Float, ?pass:Int) {
        if (ear == null) {
            return;
        }
        
        // Interlink polygon nodes in z-order
        if (pass == null && invSize != null) {
            indexCurve(ear, minX, minY, invSize);
        }
        
        var stop = ear;
        var	prev:Node = null;
        var next:Node = null;
        
        // Iterate through ears, slicing them one by one
        while (ear.prev != ear.next) {
            prev = ear.prev;
            next = ear.next;
            
            if (invSize != null ? isEarHashed(ear, minX, minY, invSize) : isEar(ear)) {
                // Cut off the triangle
                triangles.push(Std.int(prev.i / dim));
                triangles.push(Std.int(ear.i / dim));
                triangles.push(Std.int(next.i / dim));
                
                removeNode(ear);
                
                // Skipping the next vertice leads to less sliver triangles
                ear = next.next;
                stop = next.next;
                
                continue;
            }
            
            ear = next;
            
            // If we looped through the whole remaining polygon and can't find any more ears
            if (ear == stop) {
                // Try filtering points and slicing again
                if (pass == null) {
                    earcutLinked(filterPoints(ear), triangles, dim, minX, minY, invSize, 1);				
                }
                // If this didn't work, try curing all small self-intersections locally
                else if (pass == 1) {
                    ear = cureLocalIntersections(filterPoints(ear), triangles, dim);
                    earcutLinked(ear, triangles, dim, minX, minY, invSize, 2);				
                }
                // As a last resort, try splitting the remaining polygon into two
                else if (pass == 2) {
                    splitEarcut(ear, triangles, dim, minX, minY, invSize);
                }
                
                break;
            }
        }
    }

    /**
     * Check whether a polygon node forms a valid ear with adjacent nodes
     */
    static function isEar(ear:Node):Bool {
        var a = ear.prev;
        var	b = ear;
        var	c = ear.next;
        
        if (area(a, b, c) >= 0) {
            return false; // Reflex, can't be an ear
        }
        
        // Now make sure we don't have other points inside the potential ear
        var p = ear.next.next;
        
        while (p != ear.prev) {
            if (pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) && area(p.prev, p, p.next) >= 0) {
                return false;
            }
            p = p.next;
        }
        
        return true;
    }
    
    static function isEarHashed(ear:Node, minX:Float, minY:Float, invSize:Float) {
        var a = ear.prev;
        var	b = ear;
        var c = ear.next;
        
        if (area(a, b, c) >= 0) {
            return false; // Reflex, can't be an ear
        }
        
        // Triangle bbox; min & max are calculated like this for speed
        var minTX = a.x < b.x ? (a.x < c.x ? a.x : c.x) : (b.x < c.x ? b.x : c.x);
        var	minTY = a.y < b.y ? (a.y < c.y ? a.y : c.y) : (b.y < c.y ? b.y : c.y);
        var	maxTX = a.x > b.x ? (a.x > c.x ? a.x : c.x) : (b.x > c.x ? b.x : c.x);
        var	maxTY = a.y > b.y ? (a.y > c.y ? a.y : c.y) : (b.y > c.y ? b.y : c.y);
        
        // z-order range for the current triangle bbox;
        var minZ = zOrder(Std.int(minTX), Std.int(minTY), Std.int(minX), Std.int(minY), Std.int(invSize));
        var	maxZ = zOrder(Std.int(maxTX), Std.int(maxTY), Std.int(minX), Std.int(minY), Std.int(invSize));
        
        // First look for points inside the triangle in increasing z-order
        var p = ear.prevZ;
        var n = ear.nextZ;
        
        // Look for points inside the triangle in both directions
        while (p != null && p.z >= minZ && n != null && n.z <= maxZ) {
            if (p != ear.prev && p != ear.next &&
                pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) &&
                area(p.prev, p, p.next) >= 0) return false;
            p = p.prevZ;

            if (n != ear.prev && n != ear.next &&
                pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, n.x, n.y) &&
                area(n.prev, n, n.next) >= 0) return false;
            n = n.nextZ;
        }

        // Look for remaining points in decreasing z-order
        while (p != null && p.z >= minZ) {
            if (p != ear.prev && p != ear.next &&
                pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) &&
                area(p.prev, p, p.next) >= 0) return false;
            p = p.prevZ;
        }
    
        // Look for remaining points in increasing z-order
        while (n != null && n.z <= maxZ) {
            if (n != ear.prev && n != ear.next &&
                pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, n.x, n.y) &&
                area(n.prev, n, n.next) >= 0) return false;
            n = n.nextZ;
        }
    
        return true;
    }

    /**
     * Go through all polygon nodes and cure small local self-intersections
     */
    static function cureLocalIntersections(start:Node, triangles:Array<Int>, dim:Float) {
        var p = start;
        do {
            var a = p.prev;
            var	b = p.next.next;
            
            if (!equals(a, b) && intersects(a, p, p.next, b) && locallyInside(a, b) && locallyInside(b, a)) {
                triangles.push(Std.int(a.i / dim));
                triangles.push(Std.int(p.i / dim));
                triangles.push(Std.int(b.i / dim));
                
                // Remove two nodes involved
                removeNode(p);
                removeNode(p.next);
                
                p = start = b;
            }
            p = p.next;
        } 
        while (p != start);
        
        return filterPoints(p);
    }

    /**
     * Try splitting polygon into two and triangulate them independently
     */
    static function splitEarcut(start:Node, triangles:Array<Int>, dim:Float, minX:Float, minY:Float, invSize:Float) {
        // Look for a valid diagonal that divides the polygon into two
        var a = start;
        do {
            var b = a.next.next;
            while (b != a.prev) {
                if (a.i != b.i && isValidDiagonal(a, b)) {
                    // Split the polygon in two by the diagonal
                    var c = splitPolygon(a, b);
                    
                    // Filter colinear points around the cuts
                    a = filterPoints(a, a.next);
                    c = filterPoints(c, c.next);
                    
                    // Run earcut on each half
                    earcutLinked(a, triangles, dim, minX, minY, invSize);
                    earcutLinked(c, triangles, dim, minX, minY, invSize);
                    
                    return;
                }
                
                b = b.next;
            }
            
            a = a.next;
        } 
        while (a != start);
    }

    /**
     * Link every hole into the outer loop, producing a single-ring polygon without holes
     */
    static function eliminateHoles(data:Array<Float>, holeIndices:Array<Int>, outerNode:Node, dim:Int):Node {
        var queue:Array<Node> = [];
        var len:Int = holeIndices.length;
        var	start:Int = 0;
        var end:Int = 0;
        var list:Node = null;
        
        for (i in 0...len) {
            start = Std.int(holeIndices[i] * dim);
            end = i < len - 1 ? Std.int(holeIndices[i + 1] * dim) : data.length;
            list = linkedList(data, start, end, dim, false);
            
            if (list == list.next) {
                list.steiner = true;
            }
            queue.push(getLeftmost(list));
        }
        
        queue.sort(compareX);
        
        // Process holes from left to right
        for (i in 0...queue.length) {
            eliminateHole(queue[i], outerNode);
            outerNode = filterPoints(outerNode, outerNode.next);
        }
        
        return outerNode;
    }
    
    static inline function compareX(a:Node, b:Node):Int {
        return Std.int(a.x - b.x);
    }

    /**
     * Find a bridge between vertices that connects hole with an outer ring and and link it
     */
    static inline function eliminateHole(hole:Node, outerNode:Node) {
        outerNode = findHoleBridge(hole, outerNode);
        if (outerNode != null) {
            var b = splitPolygon(outerNode, hole);

            // Filter collinear points around the cuts
            filterPoints(outerNode, outerNode.next);
            filterPoints(b, b.next);
        }
    }

    /**
     * David Eberly's algorithm for finding a bridge between hole and outer polygon
     */
    static function findHoleBridge(hole:Node, outerNode:Node):Node {
        var p = outerNode;
        var hx = hole.x;
        var hy = hole.y;
        var	qx = Math.NEGATIVE_INFINITY;
        var m:Node = null;
        
        // find a segment intersected by a ray from the hole's leftmost point to the left;
        // segment's endpoint with lesser x will be potential connection point
        do {
            if (hy <= p.y && hy >= p.next.y && p.next.y != p.y) {
                var x = p.x + (hy - p.y) * (p.next.x - p.x) / (p.next.y - p.y);
                if (x <= hx && x > qx) {
                    qx = x;
                    if (x == hx) {
                        if (hy == p.y) {
                            return p;
                        }
                        if (hy == p.next.y) {
                            return p.next;
                        }
                    }
                    
                    m = p.x < p.next.x ? p : p.next;
                }
            }
            
            p = p.next;
        } 
        while (p != outerNode);
        
        if (m == null) {
            return null;
        }
        
        if (hx == qx) {
            return m.prev; // Hole touches outer segment; pick lower endpoint
        }
        
        // Look for points inside the triangle of hole point, segment intersection and endpoint;
        // If there are no points found, we have a valid connection;
        // Otherwise choose the point of the minimum angle with the ray as connection point
        
        var stop:Node = m;
        var mx:Float = m.x;
        var my:Float = m.y;
        var tanMin:Float = Math.POSITIVE_INFINITY;
        var tan:Float = 0;
        
        p = m;

        do {
            if (hx >= p.x && p.x >= mx && hx != p.x &&
                    pointInTriangle(hy < my ? hx : qx, hy, mx, my, hy < my ? qx : hx, hy, p.x, p.y)) {
    
                tan = Math.abs(hy - p.y) / (hx - p.x); // Tangential
    
                if (locallyInside(p, hole) &&
                    (tan < tanMin || (tan == tanMin && (p.x > m.x || (p.x == m.x && sectorContainsSector(m, p)))))) {
                    m = p;
                    tanMin = tan;
                }
            }
    
            p = p.next;
        } while (p != stop);
    
        return m;
    }

    /**
     * Whether sector in vertex m contains sector in vertex p in the same coordinates
     */
    static inline function sectorContainsSector(m:Node, p:Node):Bool {
        return area(m.prev, m, p.prev) < 0 && area(p.next, m, m.next) < 0;
    }

    /**
     * Interlink polygon nodes in z-order
     */
    static function indexCurve(start:Node, minX:Float, minY:Float, invSize:Float):Void {
        var p = start;
        
        do {
            if (p.z == -99999999) {
                p.z = zOrder(Std.int(p.x), Std.int(p.y), Std.int(minX), Std.int(minY), Std.int(invSize));
            }
            
            p.prevZ = p.prev;
            p.nextZ = p.next;
            p = p.next;
        } 
        while (p != start);
        
        p.prevZ.nextZ = null;
        p.prevZ = null;
        
        sortLinked(p);
    }

    /**
     * Simon Tatham's linked list merge sort algorithm
     * http://www.chiark.greenend.org.uk/~sgtatham/algorithms/listsort.html
     */
    static function sortLinked(list:Node):Node {
        var p:Node = null;
        var q:Node = null;
        var e:Node = null;
        var tail:Node;
        var numMerges:Int = 0;
        var pSize:Int = 0;
        var qSize:Int = 0;
        var	inSize = 1;
        
        do {
            p = list;
            list = null;
            tail = null;
            numMerges = 0;
            
            while (p != null) {
                numMerges++;
                q = p;
                pSize = 0;
                for (i in 0...inSize) {
                    pSize++;
                    q = q.nextZ;
                    if (q == null) {
                        break;
                    }
                }
                
                qSize = inSize;
                
                while (pSize > 0 || (qSize > 0 && q != null)) {
                    
                    if (pSize != 0 && (qSize == 0 || q == null || p.z <= q.z)) {
                        e = p;
                        p = p.nextZ;
                        pSize--;
                    }
                    else {
                        e = q;
                        q = q.nextZ;
                        qSize--;
                    }
                    
                    if (tail != null) {
                        tail.nextZ = e;
                    }
                    else {
                        list = e;
                    }
                    
                    e.prevZ = tail;
                    tail = e;
                }
                
                p = q;
            }
            
            tail.nextZ = null;
            inSize *= 2;			
        } 
        while (numMerges > 1);
        
        return list;
    }

    /**
     * z-order of a point given coords and size of the data bounding box
     */
    static function zOrder(x:Int, y:Int, minX:Int, minY:Int, invSize:Int):Int {
        // coords are transformed into non-negative 15-bit integer range
        x = Std.int(32767 * (x - minX) * invSize);
        y = Std.int(32767 * (y - minY) * invSize);
        
        x = (x | (x << 8)) & 0x00FF00FF;
        x = (x | (x << 4)) & 0x0F0F0F0F;
        x = (x | (x << 2)) & 0x33333333;
        x = (x | (x << 1)) & 0x55555555;
        
        y = (y | (y << 8)) & 0x00FF00FF;
        y = (y | (y << 4)) & 0x0F0F0F0F;
        y = (y | (y << 2)) & 0x33333333;
        y = (y | (y << 1)) & 0x55555555;
        
        return x | (y << 1);
    }

    /**
     * Find the leftmost node of a polygon ring
     */
    static function getLeftmost(start:Node):Node {
        var p = start;
        var leftmost = start;

        do {
            if (p.x < leftmost.x || (p.x == leftmost.x && p.y < leftmost.y)) leftmost = p;
            p = p.next;
        }
        while (p != start);

        return leftmost;
    }
    
    /**
     * Check if a point lies within a convex triangle
     */
    static inline function pointInTriangle(ax:Float, ay:Float, bx:Float, by:Float, cx:Float, cy:Float, px:Float, py:Float):Bool {
        return (cx - px) * (ay - py) - (ax - px) * (cy - py) >= 0 &&
               (ax - px) * (by - py) - (bx - px) * (ay - py) >= 0 &&
               (bx - px) * (cy - py) - (cx - px) * (by - py) >= 0;
    }

    /**
     * Check if a diagonal between two polygon nodes is valid (lies in polygon interior)
     */
    static inline function isValidDiagonal(a:Node, b:Node):Bool {
        return a.next.i != b.i && a.prev.i != b.i && !intersectsPolygon(a, b) && // Doesn't intersect other edges
            (locallyInside(a, b) && locallyInside(b, a) && middleInside(a, b) && // Locally visible
            (area(a.prev, a, b.prev) != 0 || area(a, b.prev, b) != 0) || // Does not create opposite-facing sectors
            equals(a, b) && area(a.prev, a, a.next) > 0 && area(b.prev, b, b.next) > 0); // Special zero-length case
    }

    /**
     * Signed area of a triangle
     */
    static inline function area(p:Node, q:Node, r:Node):Float {
        return (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    }

    /**
     * Check if two points are equal
     */
    static inline function equals(p1:Node, p2:Node):Bool {
        return return p1.x == p2.x && p1.y == p2.y;
    }

    /**
     * Check if two segments intersect
     */
    static inline function intersects(p1:Node, q1:Node, p2:Node, q2:Node):Bool {
        var o1 = sign(area(p1, q1, p2));
        var o2 = sign(area(p1, q1, q2));
        var o3 = sign(area(p2, q2, p1));
        var o4 = sign(area(p2, q2, q1));
    
        if (o1 != o2 && o3 != o4) return true; // General case
    
        if (o1 == 0 && onSegment(p1, p2, q1)) return true; // p1, q1 and p2 are collinear and p2 lies on p1q1
        if (o2 == 0 && onSegment(p1, q2, q1)) return true; // p1, q1 and q2 are collinear and q2 lies on p1q1
        if (o3 == 0 && onSegment(p2, p1, q2)) return true; // p2, q2 and p1 are collinear and p1 lies on p2q2
        if (o4 == 0 && onSegment(p2, q1, q2)) return true; // p2, q2 and q1 are collinear and q1 lies on p2q2
    
        return false;
    }

    /**
     * For collinear points p, q, r, check if point q lies on segment pr
     */
    static inline function onSegment(p:Node, q:Node, r:Node):Bool {
        return q.x <= Math.max(p.x, r.x) && q.x >= Math.min(p.x, r.x) && q.y <= Math.max(p.y, r.y) && q.y >= Math.min(p.y, r.y);
    }

    static inline function sign(num:Float):Int {
        return num > 0 ? 1 : num < 0 ? -1 : 0;
    }

    /**
     * Check if a polygon diagonal intersects any polygon segments
     */
    static function intersectsPolygon(a:Node, b:Node):Bool {
        var p = a;
        do {
            if (p.i != a.i && p.next.i != a.i && p.i != b.i && p.next.i != b.i && intersects(p, p.next, a, b)) {
                return true;
            }
            p = p.next;
        } 
        while (p != a);
        
        return false;
    }

    /**
     * Check if a polygon diagonal is locally inside the polygon
     */
    static inline function locallyInside(a:Node, b:Node) {
        return area(a.prev, a, a.next) < 0 ?
            area(a, b, a.next) >= 0 && area(a, a.prev, b) >= 0 :
            area(a, b, a.prev) < 0 || area(a, a.next, b) < 0;
    }

    /**
     * Check if the middle point of a polygon diagonal is inside the polygon
     */
    static function middleInside(a:Node, b:Node):Bool {
        var p:Node = a;
        var inside:Bool = false;
        var px:Float = (a.x + b.x) / 2;
        var py:Float = (a.y + b.y) / 2;
        
        do {
            if (((p.y > py) != (p.next.y > py)) && p.next.y != p.y &&
                    (px < (p.next.x - p.x) * (py - p.y) / (p.next.y - p.y) + p.x))
                inside = !inside;
            p = p.next;
        }
        while (p != a);
        
        return inside;
    }

    /**
     * Link two polygon vertices with a bridge; if the vertices belong to the same ring, it splits polygon into two;
     * If one belongs to the outer ring and another to a hole, it merges it into a single ring
     */
    static function splitPolygon(a:Node, b:Node):Node {
        var a2:Node = Node.get(a.i, a.x, a.y);
        var	b2:Node = Node.get(b.i, b.x, b.y);
        var	an:Node = a.next;
        var	bp:Node = b.prev;
        
        a.next = b;
        b.prev = a;
        
        a2.next = an;
        an.prev = a2;
        
        b2.next = a2;
        a2.prev = b2;
        
        bp.next = b2;
        b2.prev = bp;
        
        return b2;
    }

    /**
     * Create a node and optionally link it with previous one (in a circular doubly linked list)
     */
    static function insertNode(i:Int, x:Float, y:Float, ?last:Node):Node {
        var node = Node.get(i, x, y);
        
        if (last == null) {
            node.prev = node;
            node.next = node;
        } 
        else {
            node.next = last.next;
            node.prev = last;
            last.next.prev = node;
            last.next = node;
        }
        
        return node;
    }
    
    static function removeNode(p:Node) {
        p.next.prev = p.prev;
        p.prev.next = p.next;
        
        if (p.prevZ != null) {
            p.prevZ.nextZ = p.nextZ;
        }
        if (p.nextZ != null) {
            p.nextZ.prevZ = p.prevZ;
        }
    }
    
    static function signedArea(data:Array<Float>, start:Int, end:Int, dim:Int):Float {
        var sum:Float = 0;
        var i:Int = start;
        var j:Int = end - dim;
        while (i < end) {
            sum += (data[j] - data[i]) * (data[i + 1] + data[j + 1]);
            j = i;
            i += dim;
        }
        
        return sum;
    }
    
    // return a percentage difference between the polygon area and its triangulation area;
    // used to verify correctness of triangulation
    static function deviation(data:Array<Float>, holeIndices:Array<Int>, dim:Int, triangles:Array<Int>):Float {
        var hasHoles:Bool = holeIndices != null && holeIndices.length > 0;
        var outerLen:Int = hasHoles ? Std.int(holeIndices[0] * dim) : data.length;
        
        var polygonArea = Math.abs(signedArea(data, 0, outerLen, dim));
        if (hasHoles) {
            var len:Int = holeIndices.length;
            for (i in 0...len) {
                var start = holeIndices[i] * dim;
                var end = i < len - 1 ? holeIndices[i + 1] * dim : data.length;
                polygonArea -= Math.abs(signedArea(data, start, end, dim));
            }
        }
        
        var trianglesArea = 0.0;
        var i:Int = 0;
        while (i < triangles.length) {
            var a:Int = cast triangles[i] * dim;
            var b:Int = cast triangles[i + 1] * dim;
            var c:Int = cast triangles[i + 2] * dim;
            trianglesArea += Math.abs(
                (data[a] - data[c]) * (data[b + 1] - data[a + 1]) -
                (data[a] - data[b]) * (data[c + 1] - data[a + 1]));
                
            i += 3;
        }
        
        return polygonArea == 0 && trianglesArea == 0 ? 0 : Math.abs((trianglesArea - polygonArea) / polygonArea);
    }
    
    // turn a polygon in a multi-dimensional array form (e.g. as in GeoJSON) into a form Earcut accepts
    static public function flatten(data:Array<Array<Array<Float>>>) {
        var dim:Int = data[0][0].length;
        var vertices:Array<Float> = [];
        var holes:Array<Int> = [];
        var result:Dynamic = { };
        result.vertices = vertices;
        result.holes = holes;
        result.dimensions = dim;
        var holeIndex:Int = 0;
        
        for (i in 0...data.length) {
            for (j in 0...data[i].length) {
                for (d in 0...dim) {
                    result.vertices.push(data[i][j][d]);
                }
            }
            if (i > 0) {
                holeIndex += data[i - 1].length;
                result.holes.push(holeIndex);
            }
        }
        
        return result;
    }
    
}

class Node {
    
    public var i:Int;
    
    public var x:Float;
    public var y:Float;
    
    public var prev:Node; 
    public var next:Node;
    
    public var z:Int = -99999999;
    
    public var prevZ:Node;
    public var nextZ:Node;
    
    public var steiner:Bool;
    
    
    public function new(i:Int, x:Float, y:Float) {

        reset(i, x, y);

    }

    inline function reset(i:Int, x:Float, y:Float) {

        // vertice index in coordinates array
        this.i = i;
        
        // vertex coordinates
        this.x = x;
        this.y = y;
        
        // previous and next vertice nodes in a polygon ring
        this.prev = null;
        this.next = null;
        
        // z-order curve value
        this.z = -99999999;
        
        // previous and next nodes in z-order
        this.prevZ = null;
        this.nextZ = null;
        
        // indicates whether this is a steiner point
        this.steiner = false;

    }

/// Recycling nodes

    public static var recyclingEnabled:Bool = true;

    static var nextPoolIndex:Int = 0;

    static var pool:Array<Node> = [];

    public static function clearPool() {

        nextPoolIndex = 0;
        pool = [];

    }

    public static function recycleAllNodes() {

        nextPoolIndex = 0;

    }

    public static function get(i:Int, x:Float, y:Float) {

        if (recyclingEnabled) {
            if (nextPoolIndex == pool.length) {

                // Create new node
                var node = new Node(i, x, y);
                nextPoolIndex++;
                pool.push(node);
                return node;
            }
            else {

                // Reuse an available node
                var node = pool[nextPoolIndex];
                nextPoolIndex++;
                return node;
            }
        }
        else {
            return new Node(i, x, y);
        }

    }
    
}
