//Purpose: A file that holds the code that students fill in


//Given a ray described by an initial point P0 and a direction V both in
//world coordinates, check to see 
//if it intersects the polygon described by "vertices," an array of vec3
//values describing the location of the polygon vertices in its child frame.
//mvMatrix is a matrix describing how to transform "vertices" into world coordinates
//which you will have to do to get the correct intersection in world coordinates.
//Be sure to compute the plane normal only after you have transformed the points,
//and be sure to only compute intersections which are inside of the polygon
//(you can assume that all polygons are convex and use the area method)
function rayIntersectPolygon(P0, V, vertices, mvMatrix) {
    //TODO: Fill this in
    
    //Step 1: Make a new array of vec3s which holds "vertices" transformed 
    //to world coordinates (hint: vec3 has a function "transformMat4" which is useful)
    var transformed = [];
    for (var v = 0; v < vertices.length; v++) {
        var vtxnew = vec3.create();
        vec3.transformMat4(vtxnew, vertices[v], mvMatrix);
        transformed.push(vtxnew);
    }
    
    //Step 2: Compute the plane normal of the plane spanned by the transformed vertices
    var edge01 = vec3.create();
    var edge12 = vec3.create();
    var normal = vec3.create();
    vec3.subtract(edge01, transformed[1], transformed[0]);
    vec3.subtract(edge12, transformed[2], transformed[1]);
    vec3.cross(normal, edge01, edge12);
    
    //Step 3: Perform ray intersect plane
    var temp = vec3.create();
    vec3.subtract(temp, transformed[0], P0);
    var t = vec3.dot(temp, normal)/vec3.dot(V, normal);
    if (vec3.dot(V, normal)==0 || t <= 0) {
        //// no intersection if P0 on the face of vertices (t==0)
        return null;
    }
    var P = vec3.create();
    vec3.scaleAndAdd(P, P0, V, t);
    
    //Step 4: Check to see if the intersection point is inside of the transformed polygon
    //You can assume that the polygon is convex.  If you use the area test, you can
    //allow for some wiggle room in the two areas you're comparing (e.g. absolute difference
    //not exceeding 1e-4)
    var ref = crossProduct(transformed[transformed.length-1], transformed[0], P);
    for (var i = 0; i < transformed.length-1; i++) {
        if (vec3.dot(ref, crossProduct(transformed[i], transformed[i+1], P)) <= 0) {
            //// no intersection if P on the line of an edge (when crossProduct==0)
            return null;
        }
    } 
    
    //Step 5: Return the intersection point if it exists or null if it's outside
    //of the polygon or if the ray is perpendicular to the plane normal (no intersection)
    return {t:t, P:P}; //These are dummy values, but you should return 
    //both an intersection point and a parameter t.  The parameter t will be used to sort
    //intersections in order of occurrence to figure out which one happened first
}

//// Helper funtion to find the cross product of two vec3's AB and AC
function crossProduct(A, B, C) {
    var ab = vec3.create();
    var ac = vec3.create();
    var crs = vec3.create();
    vec3.subtract(ab, B, A);
    vec3.subtract(ac, C, A);
    vec3.cross(crs, ab, ac);
    return crs;
}


function computeBBoxes(node, mvMatrix) {
    if (node === null) {
        return;
    }
    var bbox = [];
    bbox.xmin = Infinity;
    bbox.ymin = Infinity;
    bbox.zmin = Infinity;
    bbox.xmax = -Infinity;
    bbox.ymax = -Infinity;
    bbox.zmax = -Infinity;

    if ('mesh' in node) {
        for (var f = 0; f < node.mesh.faces.length; f++) {
            var vertices = node.mesh.faces[f].getVerticesPos();
            for (var v = 0; v < vertices.length; v++) {
                var vtx = vec3.create();
                vec3.transformMat4(vtx, vertices[v], mvMatrix);
                if (vtx[0] < bbox.xmin) bbox.xmin = vtx[0];
                if (vtx[0] > bbox.xmax) bbox.xmax = vtx[0];
                if (vtx[1] < bbox.ymin) bbox.ymin = vtx[1];
                if (vtx[1] > bbox.ymax) bbox.ymax = vtx[1];
                if (vtx[2] < bbox.zmin) bbox.zmin = vtx[2];
                if (vtx[2] > bbox.zmax) bbox.zmax = vtx[2];
            }
        }
    }

    if ('children' in node) {
        for (var i = 0; i < node.children.length; i++) {
            var child = node.children[i];
            var nextmvMatrix = mat4.create();
            mat4.mul(nextmvMatrix, mvMatrix, child.transform);
            computeBBoxes(child, nextmvMatrix);

            if (child.bbox.xmin < bbox.xmin) bbox.xmin = child.bbox.xmin;
            if (child.bbox.xmax > bbox.xmax) bbox.xmax = child.bbox.xmax;
            if (child.bbox.ymin < bbox.ymin) bbox.ymin = child.bbox.ymin;
            if (child.bbox.ymax > bbox.ymax) bbox.ymax = child.bbox.ymax;
            if (child.bbox.zmin < bbox.zmin) bbox.zmin = child.bbox.zmin;
            if (child.bbox.zmax > bbox.zmax) bbox.zmax = child.bbox.zmax;
        }
    }

    //// make bounding box
    var vtx000 = vec3.fromValues(bbox.xmin, bbox.ymin, bbox.zmin);
    var vtx001 = vec3.fromValues(bbox.xmin, bbox.ymin, bbox.zmax);
    var vtx010 = vec3.fromValues(bbox.xmin, bbox.ymax, bbox.zmin);
    var vtx011 = vec3.fromValues(bbox.xmin, bbox.ymax, bbox.zmax);
    var vtx100 = vec3.fromValues(bbox.xmax, bbox.ymin, bbox.zmin);
    var vtx101 = vec3.fromValues(bbox.xmax, bbox.ymin, bbox.zmax);
    var vtx110 = vec3.fromValues(bbox.xmax, bbox.ymax, bbox.zmin);
    var vtx111 = vec3.fromValues(bbox.xmax, bbox.ymax, bbox.zmax);

    if (bbox.xmin == bbox.xmax) bbox.push([vtx000, vtx001, vtx011, vtx010]);
    if (bbox.ymin == bbox.ymax) bbox.push([vtx000, vtx001, vtx101, vtx100]);
    if (bbox.zmin == bbox.zmax) bbox.push([vtx000, vtx100, vtx110, vtx010]);
    else {
        bbox.push([vtx000, vtx001, vtx011, vtx010]);
        bbox.push([vtx000, vtx001, vtx101, vtx100]);
        bbox.push([vtx000, vtx100, vtx110, vtx010]);
        bbox.push([vtx100, vtx101, vtx111, vtx110]);
        bbox.push([vtx010, vtx011, vtx111, vtx110]);
        bbox.push([vtx001, vtx101, vtx111, vtx011]);
    }

    bbox.lines = [  [vtx000,vtx001], [vtx000,vtx010], [vtx000,vtx100], 
                    [vtx101,vtx001], [vtx101,vtx111], [vtx101,vtx100], 
                    [vtx011,vtx001], [vtx011,vtx111], [vtx011,vtx010], 
                    [vtx110,vtx111], [vtx110,vtx010], [vtx110,vtx100]];
    node.bbox = bbox;
}



function addImageSourcesFunctions(scene) {
    //Setup all of the functions that students fill in that operate directly
    //on the scene
    
    //Purpose: A recursive function provided which helps to compute intersections of rays
    //with all faces in the scene, taking into consideration the scene graph structure
    //Inputs: P0 (vec3): Ray starting point, V (vec3): ray direction
    //node (object): node in scene tree to process, 
    //mvMatrix (mat4): Matrix to put geometry in this node into world coordinates
    //excludeFace: Pointer to face object to be excluded (don't intersect with
    //the face that this point lies on)
    //Returns: null if no intersection,
    //{tmin:minimum t along ray, PMin(vec3): corresponding point, faceMin:Pointer to mesh face hit first}
    
    //NOTE: Calling this function with node = scene and an identity matrix for mvMatrix
    //will start the recursion at the top of the scene tree in world coordinates
    scene.rayIntersectFaces = function(P0, V, node, mvMatrix, excludeFace) {
        var tmin = Infinity;//The parameter along the ray of the nearest intersection
        var PMin = null;//The point of intersection corresponding to the nearest interesection
        var faceMin = null;//The face object corresponding to the nearest intersection
        if (node === null) {
            return null;
        }

        //// check bounding box
        if ('bbox' in node) {
            var intxn = null;
            for (var b = 0; b < node.bbox.length; b++) {
                intxn = rayIntersectPolygon(P0, V, node.bbox[b], mat4.create());
                if (!(intxn === null)) {
                    break;
                }
            }
            if (intxn === null) {
                return null;
            }
        }

        if ('mesh' in node) { //Make sure it's not just a dummy transformation node
            var mesh = node.mesh;
            for (var f = 0; f < mesh.faces.length; f++) {
                if (mesh.faces[f] == excludeFace) {
                    continue;//Don't re-intersect with the face this point lies on
                }
                //Intersect the ray with this polygon
                var res = rayIntersectPolygon(P0, V, mesh.faces[f].getVerticesPos(), mvMatrix);
                if (!(res === null) && (res.t < tmin)) {
                    tmin = res.t;
                    PMin = res.P;
                    faceMin = mesh.faces[f];
                }
            }
        }
        
        if ('children' in node) {
            //Recursively check the meshes of the children to make sure the ray
            //doesn't intersect any of them first
            for (var i = 0; i < node.children.length; i++) {
                var nextmvMatrix = mat4.create();
                //Multiply on the right by the next transformation of the child
                //node
                mat4.mul(nextmvMatrix, mvMatrix, node.children[i].transform);
                //Recursively intersect with the child node
                var cres = scene.rayIntersectFaces(P0, V, node.children[i], nextmvMatrix, excludeFace);
                if (!(cres === null) && (cres.tmin < tmin)) {
                    tmin = cres.tmin;
                    PMin = cres.PMin;
                    faceMin = cres.faceMin;
                }
            }
        }
        if (PMin === null) {
            return null;
        }
        return {tmin:tmin, PMin:PMin, faceMin:faceMin};
    }

    
    //Purpose: Fill in the array scene.imsources[] with a bunch of source
    //objects.  It's up to you what you put in the source objects, but at
    //the very least each object needs a field "pos" describing its position
    //in world coordinates so that the renderer knows where to draw it
    //You will certainly also need to save along pointers from an image source
    //to its parent so that when you trace paths back you know where to aim
    //Recursion is highly recommended here, since you'll be making images of 
    //images of images (etc...) reflecting across polygon faces.
    
    //Inputs: order (int) : The maximum number of bounces to take
    scene.computeImageSources = function(order) {
        console.log("debug: order is passed: " + order);
        scene.source.order = 0;//Store an order field to figure out how many 
        //bounces a particular image represents
        scene.source.rcoeff = 1.0;//Keep track of the reflection coefficient of the node that
        //gave rise to this source
        scene.source.parent = null;//Keep track of the image source's parent
        scene.source.genFace = null;//Keep track of the mesh face that generated this image
        //Remember not to reflect an image across the face that just generated it, 
        //or you'll get its parent image.  This information can also be used later
        //when tracing back paths
        scene.imsources = [scene.source];
        
        //TODO: Fill the rest of this in.  Be sure to reflect images across faces
        //in world coordinates, not the faces in the original mesh coordinates
        //See the "rayIntersectFaces" function above for an example of how to loop
        //through faces in a mesh
        if (!('polygons' in scene)) {
            scene.polygons = [];
            getFaces(scene.polygons, scene, mat4.create());
        }
        for (var n = 1; n <= order; n++) {
            var len = scene.imsources.length;
            for (var k = 0; k < len; k++) {
                var source = scene.imsources[k];
                for (var i = 0; i < scene.polygons.length; i++) {
                    var f = scene.polygons[i];
                    if (f.face != source.genFace) {
                        var normal = f.face.getNormal();
                        var trans = mat3.create();
                        mat3.normalFromMat4(trans, f.matrix);
                        vec3.transformMat3(normal, normal, trans);
                        f.face.normal = normal;//// 
                        var vtx0 = vec3.create();
                        vec3.transformMat4(vtx0, f.face.getVerticesPos()[0], f.matrix);
                        var temp = vec3.create();
                        vec3.subtract(temp, vtx0, source.pos);
                        var t = vec3.dot(temp, normal)/vec3.dot(normal, normal);
                        var P = vec3.create();
                        vec3.scaleAndAdd(P, source.pos, normal, 2*t);//// what if image on a face (t==0)?
                        image = {pos:P};
                        image.order = n;
                        image.rcoeff = source.rcoeff * f.rcoeff;
                        image.parent = source;
                        image.genFace = f.face;
                        scene.imsources.push(image);
                    }
                }
            }
        }
    }

    function getFaces(polygons, node, mvMatrix) {
        if (node === null) {
            return;
        }
        if ('mesh' in node) {
            var mesh = node.mesh;
            for (var f = 0; f < mesh.faces.length; f++) {
                polygons.push({face:mesh.faces[f], rcoeff:node.rcoeff, matrix:mvMatrix});
            }
        }
        if ('children' in node) {
            for (var i = 0; i < node.children.length; i++) {
                var nextmvMatrix = mat4.create();
                mat4.mul(nextmvMatrix, mvMatrix, node.children[i].transform);
                getFaces(polygons, node.children[i], nextmvMatrix);
            }
        }
    }
    
    scene.extractPathsBinaural = function() {
        console.log("extracting binaural paths");
        var enableBBox = true;//// change to false for testing execution time
        var begin = (new Date()).getTime();
        var text = "s to extract paths without using bounding boxes";
        if (enableBBox) {           
            text = "s to extract paths with existing bounding boxes";
            if (!('bbox' in scene)) {
                text = "s to compute bounding boxes and extract paths";
                computeBBoxes(scene, mat4.create());
                //// Log: strangely, I can't compute bboxes right after adding scene funcrions.
                //// mesh.faces is empty after parseNode() and even after everything's set up.
            }
        }
        scene.pathsL = [];
        scene.pathsR = [];
        var leftEar = vec3.create();
        var rightEar = vec3.create();
        var earWidth = 5;
        vec3.scaleAndAdd(leftEar,scene.receiver.pos,scene.receiver.right,-earWidth);
        vec3.scaleAndAdd(rightEar,scene.receiver.pos,scene.receiver.right,earWidth);

        //left
        loop:
        for (var i = 0; i < scene.imsources.length; i++) {
            var path = [{pos:leftEar,dummy:null}];//debug
            var src = scene.imsources[i];
            var order = src.order;
            var base = leftEar; //debug
            var ray = vec3.create();
            vec3.subtract(ray, src.pos, base);
            var excludeFace = null;
            for (var g = 0; g < order; g++) {
                var intxn = scene.rayIntersectFaces(base, ray, scene, mat4.create(), excludeFace);
                if (intxn === null || intxn.faceMin != src.genFace || intxn.tmin >= 1) {
                    //// continue if image on its genFace (tmin == 1)
                    // transFlag = 1 when transmission is on
                    continue loop;   
                }
                base = intxn.PMin;

                //// debugging feature
                var normal = vec3.clone(src.genFace.normal);
                vec3.scale(normal, normal, vec3.dot(normal, ray));
                vec3.normalize(normal, normal);
                vec3.scale(normal, normal, 0.3);
                var head = vec3.create();
                vec3.subtract(head, base, normal);

                path.push({pos:base, rcoeff:src.rcoeff, dir:head});
                excludeFace = src.genFace;
                src = src.parent;
                vec3.subtract(ray, src.pos, base);
            }
            intxn = scene.rayIntersectFaces(base, ray, scene, mat4.create(), excludeFace);
            if (intxn === null || intxn.tmin >= 1) {
                //// OK if source on a face (tmin == 1)
                path.push(src);
                scene.pathsL.push(path);
            }
        }

        //right
        loop:
        for (var i = 0; i < scene.imsources.length; i++) {
            var path = [{pos:rightEar, dummy:null}];//debug
            var src = scene.imsources[i];
            var order = src.order;
            var base = rightEar; //debug
            var ray = vec3.create();
            vec3.subtract(ray, src.pos, base);
            var excludeFace = null;
            for (var g = 0; g < order; g++) {
                var intxn = scene.rayIntersectFaces(base, ray, scene, mat4.create(), excludeFace);
                if (intxn === null || intxn.faceMin != src.genFace || intxn.tmin >= 1) {
                    //// continue if image on its genFace (tmin == 1)
                    // transFlag = 1 when transmission is on
                    continue loop;   
                }
                base = intxn.PMin;

                //// debugging feature
                var normal = vec3.clone(src.genFace.normal);
                vec3.scale(normal, normal, vec3.dot(normal, ray));
                vec3.normalize(normal, normal);
                vec3.scale(normal, normal, 0.3);
                var head = vec3.create();
                vec3.subtract(head, base, normal);

                path.push({pos:base, rcoeff:src.rcoeff, dir:head});
                excludeFace = src.genFace;
                src = src.parent;
                vec3.subtract(ray, src.pos, base);
            }
            intxn = scene.rayIntersectFaces(base, ray, scene, mat4.create(), excludeFace);
            if (intxn === null || intxn.tmin >= 1) {
                //// OK if source on a face (tmin == 1)
                path.push(src);
                scene.pathsR.push(path);
            }
        }

        var end = (new Date()).getTime();
        console.log((end-begin)/1000.0 + text);
    }

    //Purpose: Based on the extracted image sources, trace back paths from the
    //receiver to the source, checking to make sure there are no occlusions
    //along the way.  Remember, you're always starting by tracing a path from
    //the receiver to the image, and then from the intersection point with
    //that image's corresponding face to the image's parent, and so on
    //all the way until you get back to the original source.
    
    //Fill in the array scene.paths, where each element of the array is itself
    //an array of objects describing vertices along the path, starting
    //with the receiver and ending with the source.  Each object in each path
    //array should contain a field "pos" which describes the position, as well
    //as an element "rcoeff" which stores the reflection coefficient at that
    //part of the path, which will be used to compute decays in "computeInpulseResponse()"
    //Don't forget the direct path from source to receiver!
    scene.extractPaths = function(flag) {
        console.log("extracting path with flag: " + flag);
        var enableBBox = true;//// change to false for testing execution time
        var begin = (new Date()).getTime();
        var text = "s to extract paths without using bounding boxes";
        if (enableBBox) {           
            text = "s to extract paths with existing bounding boxes";
            if (!('bbox' in scene)) {
                text = "s to compute bounding boxes and extract paths";
                computeBBoxes(scene, mat4.create());
                //// Log: strangely, I can't compute bboxes right after adding scene funcrions.
                //// mesh.faces is empty after parseNode() and even after everything's set up.
            }
        }
        scene.paths = [];
        
        //TODO: Finish this. Extract the rest of the paths by backtracing from
        //the image sources you calculated.  Return an array of arrays in
        //scene.paths.  Recursion is highly recommended
        //Each path should start at the receiver and end at the source
        //(or vice versa), so scene.receiver should be the first element 
        //and scene.source should be the last element of every array in 
        //scene.paths
        loop:
        for (var i = 0; i < scene.imsources.length; i++) {
            var path = [scene.receiver];
            var src = scene.imsources[i];
            var order = src.order;
            var base = scene.receiver.pos;
            var ray = vec3.create();
            vec3.subtract(ray, src.pos, base);
            var excludeFace = null;
            for (var g = 0; g < order; g++) {
                var intxn = scene.rayIntersectFaces(base, ray, scene, mat4.create(), excludeFace);
                if (intxn === null || intxn.faceMin != src.genFace || intxn.tmin >= 1) {
                    //// continue if image on its genFace (tmin == 1)
                    // transFlag = 1 when transmission is on
                    if (flag==0) continue loop;   
                }
                base = intxn.PMin;

                //// debugging feature
                var normal = vec3.clone(src.genFace.normal);
                vec3.scale(normal, normal, vec3.dot(normal, ray));
                vec3.normalize(normal, normal);
                vec3.scale(normal, normal, 0.3);
                var head = vec3.create();
                vec3.subtract(head, base, normal);

                path.push({pos:base, rcoeff:src.rcoeff, dir:head});
                excludeFace = src.genFace;
                src = src.parent;
                vec3.subtract(ray, src.pos, base);
            }
            intxn = scene.rayIntersectFaces(base, ray, scene, mat4.create(), excludeFace);
            if (intxn === null || intxn.tmin >= 1) {
                //// OK if source on a face (tmin == 1)
                path.push(src);
                scene.paths.push(path);
            }
        }
        var end = (new Date()).getTime();
        console.log((end-begin)/1000.0 + text);
    }
    
    //Binaural sound
    scene.computeImpulseResponseBinaural = function(Fs) {
        var SVel = 340;//Sound travels at 340 meters/second

        scene.impulsesL = [];
        scene.impulsesR = [];
        var time_maxL = 0;
        var time_maxR = 0;
        var p = 0.8;
        //left
        loop:
        for (var f = 0; f < scene.pathsL.length; f++) { // for each path
            var path = scene.pathsL[f];
            var rec = path[0];
            var src = path[path.length-1];
            var locDis = Math.sqrt(vec3.squaredDistance(rec.pos, path[1].pos));
            var totalDis = locDis;
            var atten = 1;
            atten *= attenuate(p,locDis);
            for (var j = 1; j < path.length-1; j++) { // for each bounces, specially handle first and last segments
                if (j==path.length-2){
                    locDis = Math.sqrt(vec3.squaredDistance(path[j].pos,src.pos));
                }
                else locDis = Math.sqrt(vec3.squaredDistance(path[j].pos,path[j+1].pos));
                totalDis+=locDis;
                atten *= path[j].rcoeff*attenuate(p,locDis);
            }
            var time = totalDis/SVel;
            // console.log("Time is: "+time);
            scene.impulsesL.push({time:time,atten:atten})
            if (time>time_maxL) time_maxL=time;

        }
        
        scene.impulseRespL=new Float32Array(Math.ceil(time_maxL*Fs)); 
        for (var i=0;i<scene.impulsesL.length;i++){
            var ind = findNear(scene.impulsesL[i].time*Fs);
            scene.impulseRespL[ind] += scene.impulsesL[i].atten;    
        }

        //right
        loop:
        for (var f = 0; f < scene.pathsR.length; f++) { // for each path
            var path = scene.pathsR[f];
            var rec = path[0];
            var src = path[path.length-1];
            var locDis = Math.sqrt(vec3.squaredDistance(rec.pos, path[1].pos));
            var totalDis = locDis;
            var atten = 1;
            atten *= attenuate(p,locDis);
            for (var j = 1; j < path.length-1; j++) { // for each bounces, specially handle first and last segments
                if (j==path.length-2){
                    locDis = Math.sqrt(vec3.squaredDistance(path[j].pos,src.pos));
                }
                else locDis = Math.sqrt(vec3.squaredDistance(path[j].pos,path[j+1].pos));
                totalDis+=locDis;
                atten *= path[j].rcoeff*attenuate(p,locDis);
            }
            var time = totalDis/SVel;
            // console.log("Time is: "+time);
            scene.impulsesR.push({time:time,atten:atten})
            if (time>time_maxR) time_maxR=time;

        }
        
        scene.impulseRespR=new Float32Array(Math.ceil(time_maxR*Fs)); 
        for (var i=0;i<scene.impulsesR.length;i++){
            var ind = findNear(scene.impulsesR[i].time*Fs);
            scene.impulseRespR[ind] += scene.impulsesR[i].atten;    
        }
    }

    //Inputs: Fs: Sampling rate (samples per second)
    scene.computeImpulseResponse = function(Fs) {
        var SVel = 340;//Sound travels at 340 meters/second
        //TODO: Finish this.  Be sure to scale each bounce by 1/(1+r^p), 
        //where r is the length of the line segment of that bounce in meters
        //and p is some integer less than 1 (make it smaller if you want the 
        //paths to attenuate less and to be more echo-y as they propagate)
        //Also be sure to scale by the reflection coefficient of each material
        //bounce (you should have stored this in extractPaths() if you followed
        //those directions).  Use some form of interpolation to spread an impulse
        //which doesn't fall directly in a bin to nearby bins
        //Save the result into the array scene.impulseResp[]

        scene.impulses = [];
        var time_max = 0;
        loop:
        for (var f = 0; f < scene.paths.length; f++) { // for each path
            var p = 0.8;
            var path = scene.paths[f];
            var rec = path[0];
            var src = path[path.length-1];
            var locDis = Math.sqrt(vec3.squaredDistance(rec.pos, path[1].pos));
            var totalDis = locDis;
            var atten = 1;
            atten *= attenuate(p,locDis);
            for (var j = 1; j < path.length-1; j++) { // for each bounces, specially handle first and last segments
                if (j==path.length-2){
                    locDis = Math.sqrt(vec3.squaredDistance(path[j].pos,src.pos));
                }
                else locDis = Math.sqrt(vec3.squaredDistance(path[j].pos,path[j+1].pos));
                totalDis+=locDis;
                atten *= path[j].rcoeff*attenuate(p,locDis);
            }
            var time = totalDis/SVel;
            // console.log("Time is: "+time);
            scene.impulses.push({time:time,atten:atten})
            if (time>time_max) time_max=time;

        }
        
        scene.impulseResp=new Float32Array(Math.ceil(time_max*Fs)); 
        for (var i=0;i<scene.impulses.length;i++){
            var ind = findNear(scene.impulses[i].time*Fs);
            scene.impulseResp[ind] += scene.impulses[i].atten;    
        }
    }

    //helper calculate attenuation
    function attenuate(p,dis){
        return 1/((1+dis)^p);
    }
    //helper find nearest bin
    function findNear(num){
        var ceil = Math.ceil(num);
        if ((ceil-num)<=0.5) return ceil;
        else return Math.floor(num);
    }
}


function buildCityA(scene) {
    var m = 1; //scale factor
    var d = 500; //city floor dimension
    var smin = 45; //street width min
    var smax = 75; //street width max
    var bmin = 30; //building size min
    var bmax = 50; //building size max
    var hmin = 100; //building height min
    var hmax = 250; //building height max

    var bx = [];
    var sx = [];
    var dx = d;
    while (dx > 0) {
        var b = Math.floor((Math.random() * (bmax-bmin+1)) + bmin);
        var s = Math.floor((Math.random() * (smax-smin+1)) + smin);
        if (dx < b+s) {
            if (dx < bmin) {
                break;
            }
            else if (dx < b) {
                b = dx;
            }
            s = dx-b;
        }
        bx.push(b);
        sx.push(s); 
        dx -= b+s;
    }

    var by = [];
    var sy = [];
    var dy = d;
    while (dy > 0) {
        var b = Math.floor((Math.random() * (bmax-bmin+1)) + bmin);
        var s = Math.floor((Math.random() * (smax-smin+1)) + smin);
        if (dy < b+s) {
            if (dy < bmin) {
                break;
            }
            else if (dy < b) {
                b = dy;
            }
            s = dy-b;
        }
        by.push(b);
        sy.push(s); 
        dy -= b+s;
    }

    var children = [{mesh:"meshes/square.off",
                    color:[1, 1, 1],
                    rcoeff:0.5,
                    transform:[ d, 0, 0, 0,
                                0, 1, 0, 0,
                                0, 0, d, 0,
                                0, 0, 0, 1]}];
    var x = -d/2;
    for (var i = 0; i < bx.length; i++) {
        var y = -d/2;
        x += bx[i];
        for (var j = 0; j < by.length; j++) {
            y += by[j];
            var h = Math.floor((Math.random() * (hmax-hmin+1)) + hmin);
            console.log(bx[i] + " by " + by[j] + " building centered at (" + (x-bx[i]/2) + ", " + (y-by[j]/2) + ")");
            children.push({ mesh:"meshes/box.off",
                            color:[Math.random(), Math.random(), Math.random()],
                            rcoeff:0.5,
                            transform:[ bx[i], 0, 0, x-bx[i]/2,
                                        0, h, 0, h/2,
                                        0, 0, by[j], y-by[j]/2,
                                        0, 0, 0, 1]});
            y += sy[j];
        }
        x += sx[i];
    }

    var rcvx = -d/2;
    var rcvy = -d/2;
    var rx = Math.floor(Math.random() * (bx.length-1));
    for (var k = 0; k <= rx; k++) {
        rcvx += bx[k] + sx[k];
    }
    var srcx = rcvx + bx[rx+1] + sx[rx+1]/2;//// 
    rcvx -= sx[rx]/2;

    var ry = Math.floor(Math.random() * (by.length-1));
    for (var k = 0; k <= ry; k++) {
        rcvy += by[k] + sy[k];
    }
    var srcy = rcvy + by[ry+1] + sy[ry+1]/2;//// 
    rcvy -= sy[ry]/2;
    if (Math.random() < 0.5) {
        srcy = rcvy;
    }

    // var srcx = rcvx + bx[rx+1] + sx[rx+1]/2;
    // var srcy = -d/2;
    // var lx = Math.floor(Math.random() * bx.length);
    // for (var k = 0; k <= lx; k++) {
    //     srcx += bx[k] + sx[k];
    // }
    // srcx -= sx[lx]/2;
    // var ly = Math.floor(Math.random() * by.length);
    // for (var k = 0; k <= ly; k++) {
    //     srcy += by[k] + sy[k];
    // }
    // srcy -= sy[ly]/2;

    var parent = {children:children, transform:[m, 0, 0, 0,
                                                0, m, 0, -hmax*m,
                                                0, 0, m, 0,
                                                0, 0, 0, 1]};
    return {children: [parent], receiver: [rcvx*m, 1.5-hmax*m, rcvy*m], source: [srcx*m, 1.5-hmax*m, srcy*m]};
}
