.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


polymake and Optimization
=========================

By Sebastian Pokutta and Marc E. Pfetsch.

Introduction
------------

Polymake offers many interesting features that can help (discrete)
optimizers to analyze optimization problems. For example

-  linear optimization can be performed exactly and visualized in small
   dimensions

-  the convex hull of feasible points of an integer program can be
   computed and analyzed

-  Hilbert bases can be computed

There are several other tutorials that cover similar topics:

-  `ILP and Hilbert Bases <ilp_and_hilbertbases>`__

-  `Gomory Cuts <empty/michaels_tutorial2>`__

-  `Branch and Bound <empty/michaels_tutorial>`__

-  `Tutorial for Lattice Polytopes <lattice_polytopes_tutorial>`__

This tutorial is targeted towards the optimization community, since,
surprisingly, polymake does not seem to be well known in this community.
In particular, the community still tends to use the quite old program
``porta`` to compute convex hulls or inequality descriptions. While
``porta`` still does a decent job here, ``polymake`` offers a much
broader feature set. Polymake supports several convex hull algorithms
which might be better suited depending on the data. Moreover it offers
many visualization tools that can help to better *understand* a given
polytope. We think that polymake has many advantages for discrete
optimizers and hope that this tutorial will help to spread the usage of
polymake.

You can find files of the example LPs in the folder
demo/files/optimization in your polymake directory.

Input: lp2poly
--------------

The first important step is to get the desired input into ``polymake``.
While it is possible to define polytopes in ``polymake`` right away, a
more convenient methods, especially for people working in discrete
optimization, is the possibility to read an optimization problem from a
file. This has the advantage that input that is usually generated
outside ``polymake`` can simply be loaded. A standard file format for
linear (integer) programs is the so-called CPLEX LP-format, which is
well known in the optimization community and can be read into
``polymake`` in a very straightforward way as shown in the following
example:

::

   Minimize
    obj:  x1 + x2 + x3
   Subject to
    C1: x1 + x2 + x3 <= 2
   Bounds
    0 <= x1 <= 1
    0 <= x2 <= 1
    0 <= x3 <= 1
   End

Thus, the file describes a 0/1-cube in three dimensions. It should be
easy to adapt this format to other cases (If for example ``x1`` does not
have any bounds you can write ``x1 free`` instead).

Now assume that this example is contained in file ``c3t.lp``. We create
a polytope from the file via:


::

    polymake> $f=lp2poly('files/optimization/c3t.lp');

The polytope ``$f`` is coded via floating point numbers:


::

    polymake> print $f->type->full_name;
    Polytope<Rational>




We convert it to a rational polytope via:


::

    polymake> $p = new Polytope<Rational>($f);

Now, ``$p`` contains data for the linear optimization problem in the
file ``ct3.lp``. The corresponding polyhedron is stored directly in
``$p``. The information about the linear objective is contained in
``$p->LP``.

Linear Optimization
-------------------

Polymake can be used to perform several actions related to linear
optimization (linear programming - LP). For instance, one can exactly
solve a linear program (via lrs or cdd). Before we explain the
corresponding usage, we first need to have a linear optimization problem
at hand.

Assuming that we are given the above example in variable ``$p``, the
linear objective can be accessed via


::

    polymake> print $p->LP->LINEAR_OBJECTIVE;
    0 1 1 1




Thus - as described in the file - the objective function coefficients
are 1 for all three variables (and there is an offset of 0).

Now, we can solve the corresponding linear program via


::

    polymake> print $p->LP->MAXIMAL_VALUE;
    2




Thus, the maximal value that we can obtain via the above linear
objective function is 2. We can also get an optimal vertex via


::

    polymake> print $p->LP->MAXIMAL_VERTEX;
    1 0 1 1




This vertex corresponds to setting ``x1=0, x2=1, x3=1``. The optimal
face can also be computed:


::

    polymake> print $p->LP->MAXIMAL_FACE;
    {4 5 6}




This means that the optimal face is the convex hull of three vertices
(with indices 4, 5, 6).

Of course, by replacing ``MAXIMAL`` by ``MINIMAL``, we obtain the
corresponding values for minimization.

The directed graph obtained by directing the graph of the polytope in
the direction of increasing objective function can be obtained via


::

    polymake> $p->VISUAL->DIRECTED_GRAPH;


.. raw:: html

    <!--
    polymake for andrew
    Thu Mar 28 11:14:51 2019
    c3t
    -->
    
    
    <html>
       <head>
          <title>c3t</title>
          <style>
    /*
    // COMMON_CODE_BLOCK_BEGIN
    */
             html{overflow: scroll;}
             body { font-family: Arial, Helvetica, sans-serif}
             strong{font-size: 18px;}
             canvas { z-index: 8; }
             input[type='range'] {}
             input[type='radio'] {margin-left:0;}
             input[type='checkbox'] {margin-right:7px; margin-left: 0px; padding-left:0px;}
             .group{padding-bottom: 40px;}
             .settings * {z-index: 11; }
             .settings{z-index: 10; margin-left: 30px; display: none; width: 14em; height: 90%; border: solid 1px silver; padding: 2px; overflow-y: scroll; background-color: white }
             .indented{margin-left: 20px; margin-top: 15px; padding-bottom: 0px;} 
             .shownObjectsList{overflow: auto; max-width: 150px; max-height: 150px;}
             .showSettingsButton{display: block; z-index: 12; position: absolute }
             .hideSettingsButton{display: none; z-index: 12; position: absolute; opacity: 0.5}
             .resetButton{margin-top: 20px;}
             button{margin-left: 0;}
             img{cursor: pointer;}
             .suboption{padding-top: 30px;}
             .transparency{display: none;}
             .labelsCheckbox{margin-top: 10px;}
    
    
             input[type=range] {
               -webkit-appearance: none;
               padding:0; 
               width:90%; 
               margin-left: auto;
               margin-right: auto;
               margin-top: 20px;
               display: block;	
             }
             input[type=range]:focus {
               outline: none;
             }
             input[type=range]::-webkit-slider-runnable-track {
               height: 4px;
               cursor: pointer;
               animate: 0.2s;
               box-shadow: 0px 0px 0px #000000;
               background: #E3E3E3;
               border-radius: 0px;
               border: 0px solid #000000;
             }
             input[type=range]::-webkit-slider-thumb {
               box-shadow: 1px 1px 2px #B8B8B8;
               border: 1px solid #ABABAB;
               height: 13px;
               width: 25px;
               border-radius: 20px;
               background: #E0E0E0;
               cursor: pointer;
               -webkit-appearance: none;
               margin-top: -5px;
             }
             input[type=range]:focus::-webkit-slider-runnable-track {
               background: #E3E3E3;
             }
             input[type=range]::-moz-range-track {
               height: 4px;
               cursor: pointer;
               animate: 0.2s;
               box-shadow: 0px 0px 0px #000000;
               background: #E3E3E3;
               border-radius: 0px;
               border: 0px solid #000000;
             }
             input[type=range]::-moz-range-thumb {
               box-shadow: 1px 1px 2px #B8B8B8;
               border: 1px solid #ABABAB;
               height: 13px;
               width: 25px;
               border-radius: 20px;
               background: #E0E0E0;
               cursor: pointer;
             }
             input[type=range]::-ms-track {
               height: 4px;
               cursor: pointer;
               animate: 0.2s;
               background: transparent;
               border-color: transparent;
               color: transparent;
             }
             input[type=range]::-ms-fill-lower {
               background: #E3E3E3;
               border: 0px solid #000000;
               border-radius: 0px;
               box-shadow: 0px 0px 0px #000000;
             }
             input[type=range]::-ms-fill-upper {
               background: #E3E3E3;
               border: 0px solid #000000;
               border-radius: 0px;
               box-shadow: 0px 0px 0px #000000;
             }
             input[type=range]::-ms-thumb {
               box-shadow: 1px 1px 2px #B8B8B8;
               border: 1px solid #ABABAB;
               height: 13px;
               width: 25px;
               border-radius: 20px;
               background: #E0E0E0;
               cursor: pointer;
             }
             input[type=range]:focus::-ms-fill-lower {
               background: #E3E3E3;
             }
             input[type=range]:focus::-ms-fill-upper {
               background: #E3E3E3;
             }
    /*
    // COMMON_CODE_BLOCK_END
    */
    		</style>
       </head>
    
    <body>
    
    		<div id='settings_6' class='settings'>
    			<div class=group id='explode_6'>
    				<strong>Explode</strong>
    				<input id='explodeRange_6' type='range' min=0 max=6 step=0.01 value=0>
    				<div class=indented><input id='explodeCheckbox_6' type='checkbox'>Automatic explosion</div>
    				<div class=suboption>Exploding speed</div>
    				<input id='explodingSpeedRange_6' type='range' min=0 max=0.5 step=0.001 value=0.05>
    			</div>
    
    			
    			<div class=group id='transparency_6' class='transparency'>
    				<strong>Transparency</strong>
    				<input id='transparencyRange_6' type='range' min=0 max=1 step=0.01 value=0>
    			</div>
    			
    			<div class=group id='rotation_6'>
    				<strong>Rotation</strong>
    				<div class=indented>
    					<div><input type='checkbox' id='changeRotationX_6'> x-axis</div>
    					<div><input type='checkbox' id='changeRotationY_6'> y-axis</div>
    					<div><input type='checkbox' id='changeRotationZ_6'> z-axis</div>
    					<button id='resetButton_6' class='resetButton' >Reset</button>
    				</div>
    
    				<div class=suboption>Rotation speed</div>
    				<input id='rotationSpeedRange_6' type='range' min=0 max=5 step=0.01 value=2>
    
    			</div>
    
    
    			<div class=group id='display_6'>
    				<strong>Display</strong>
    				<div class=indented>
    					<div id='shownObjectsList_6' class='shownObjectsList'></div>
    					<div class='labelsCheckbox'><input type='checkbox' id='labelsCheckboxInput_6' checked>Labels</div>
    				</div>
    			</div>
    
    
    			<div class=group id='svg_6'>
    				<strong>SVG</strong>
    				<div class=indented>
    					<form>
    						<input type="radio" name='screenshotMode' value='download' id='download_6' checked> Download<br>
    						<input type="radio" name='screenshotMode' value='tab' id='tab_6' > New tab<br>
    					</form>
    					<button id='takeScreenshot_6'>Screenshot</button>
    				</div>
    			</div>
    
    		</div>	<!-- end of settings -->
    		<img id='hideSettingsButton_6' style="display: none" class='hideSettingsButton' src='/kernelspecs/polymake/close.svg' width=20px">
    		<img id='showSettingsButton_6' class='showSettingsButton' src='/kernelspecs/polymake/menu.svg' width=20px">
    <div id="model14578961617"></div>
    
    <script>
    requirejs.config({
      paths: {
        three: '/kernelspecs/polymake/three',
        Detector: '/kernelspecs/polymake/Detector',
        SVGRenderer: '/kernelspecs/polymake/SVGRenderer',
        CanvasRenderer: '/kernelspecs/polymake/CanvasRenderer',
        Projector: '/kernelspecs/polymake/Projector',
        TrackballControls: '/kernelspecs/polymake/TrackballControls'
      },
      shim: {
        'three':
        {
          exports: 'THREE'
        },
        'Detector':
        {
          deps: [ 'three' ],
          exports: 'Detector'
        },
        'SVGRenderer':
        {
          deps: [ 'three' ],
          exports: 'THREE.SVGRenderer'
        },
        'CanvasRenderer':
        {
          deps: [ 'three' ],
          exports: 'THREE.CanvasRenderer'
        },
        'Projector':
        {
          deps: [ 'three' ],
          exports: 'THREE.Projector'
        },
        'TrackballControls':
        {
          deps: [ 'three' ],
          exports: 'THREE.TrackballControls'
        }
      }
    });
    require(['three'],function(THREE){
        window.THREE = THREE;
      require(['Detector','SVGRenderer','CanvasRenderer','Projector','TrackballControls'],function(Detector,SVGRenderer,CanvasRenderer,Projector,TrackballControls){
          THREE.SVGRenderer = SVGRenderer;
          THREE.CanvasRenderer = CanvasRenderer;
          THREE.Projector = Projector;
          THREE.TrackballControls = TrackballControls;
    
    // COMMON_CODE_BLOCK_BEGIN
    	var foldable = false;
       var container = document.getElementById( 'model14578961617' );
       var renderer = Detector.webgl? new THREE.WebGLRenderer({antialias: true}): new THREE.CanvasRenderer({antialias: true});
    	var svgRenderer = new THREE.SVGRenderer({antialias: true});
                var box = document.getElementsByClassName( 'output_subarea' )[0];
             var notebook = document.getElementById( 'notebook_panel' );
    
       var width = box.clientWidth - 25;
       var height = notebook.clientHeight * 0.8;
       renderer.setSize(width, height);
       svgRenderer.setSize(width, height);
       renderer.setClearColor(0xFFFFFF, 1);
       svgRenderer.setClearColor(0xFFFFFF, 1);
    
       container.appendChild(renderer.domElement);
    
       var scene = new THREE.Scene();
       var camera = new THREE.PerspectiveCamera(75, width/height, 0.1, 1000);
    
       var renderid;
    
       camera.position.set(0, 0, 5);
       camera.lookAt(0, 0, 0);
       camera.up.set(0, 1, 0);
    
       // class to allow move points together with labels and spheres
       var PMPoint = function (x,y,z) {
          this.vector = new THREE.Vector3(x,y,z);
          this.sprite = null;
          this.sphere = null;
       }
       PMPoint.prototype.makelabel = function(label) {
          this.sprite = textSprite( label );
          this.sprite.position.copy(this.vector);
       }
       PMPoint.prototype.makesphere = function(radius,material) {
          this.sphere = new THREE.Mesh(new THREE.SphereGeometry(radius), material);
          this.sphere.position.copy(this.vector);
       }
    
       PMPoint.prototype.setX = function(x) {
          this.vector.setX(x);
          if (this.sprite) {
             this.sprite.position.setX(x);
          }
          if (this.sphere) {
             this.sphere.position.setX(x);
          }
       };
       PMPoint.prototype.setY = function(y) {
          this.vector.setY(y);
          if (this.sprite) {
             this.sprite.position.setY(y);
          }
          if (this.sphere) {
             this.sphere.position.setY(y);
          }
       };
       PMPoint.prototype.setZ = function(z) {
          this.vector.setZ(z);
          if (this.sprite) {
             this.sprite.position.setZ(z);
          }
          if (this.sphere) {
             this.sphere.position.setZ(z);
          }
       };
       PMPoint.prototype.set = function(x,y,z) {
          this.vector.set(x,y,z);
          if (this.sprite) {
             this.sprite.position.set(x,y,z);
          }
          if (this.sphere) {
             this.sphere.position.set(x,y,z);
          }
       };
       PMPoint.prototype.add = function(o) {
          if (this.sprite) {
             o.add(this.sprite);
          }
          if (this.sphere) {
             o.add(this.sphere);
          }
       };
    
    
       var controls = new THREE.TrackballControls(camera, container);
    	controls.zoomSpeed = 0.2;
    	controls.rotateSpeed = 4;
    
       var all_objects = [];
       var centroids = [];
       // select the target node
       var target = document.querySelector('#model14578961617');
    
       // create an observer instance
       var observer = new MutationObserver(function(mutations) {
          mutations.forEach(function(mutation) {
             if (mutation.removedNodes && mutation.removedNodes.length > 0) {
                cancelAnimationFrame(renderId);
                observer.disconnect();
                console.log("cancelled frame "+renderId);
             }
          });
       });
    
       // configuration of the observer:
       var config = { childList: true, characterData: true }
    
       // pass in the target node, as well as the observer options
       while (target) {
          if (target.className=="output") {
             observer.observe(target, config);
             break;
          }
          target = target.parentNode;
       }
    
    // COMMON_CODE_BLOCK_END
    
       var objectnames = ["c3t","GRAPH directed with LP unnamed#0"];
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 0, 1));
       allpoints.push(new PMPoint(1, 0, 0));
       allpoints.push(new PMPoint(0, 0, 0));
       allpoints.push(new PMPoint(0, 1, 0));
       allpoints.push(new PMPoint(1, 1, 0));
       allpoints.push(new PMPoint(1, 0, 1));
       allpoints.push(new PMPoint(0, 1, 1));
    
       <!-- Vertex style -->
       var points_material = new THREE.MeshBasicMaterial ( {color: 0xFF0000, } );
    
       points_material.side = THREE.DoubleSide;
       points_material.transparent = true;
    
       <!-- POINTS -->
       allpoints[0].makesphere(0.02,points_material);
       allpoints[1].makesphere(0.02,points_material);
       allpoints[2].makesphere(0.02,points_material);
       allpoints[3].makesphere(0.02,points_material);
       allpoints[4].makesphere(0.02,points_material);
       allpoints[5].makesphere(0.02,points_material);
       allpoints[6].makesphere(0.02,points_material);
       allpoints[0].makelabel("0");
       allpoints[1].makelabel("1");
       allpoints[2].makelabel("2");
       allpoints[3].makelabel("3");
       allpoints[4].makelabel("4");
       allpoints[5].makelabel("5");
       allpoints[6].makelabel("6");
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
       faces.vertices.push(allpoints[3].vector);
       faces.vertices.push(allpoints[4].vector);
       faces.vertices.push(allpoints[5].vector);
       faces.vertices.push(allpoints[6].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(5, 4, 6, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(6, 3, 2, undefined, undefined, 0));
       faces.faces.push(new THREE.Face3(6, 2, 0, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(1, 4, 5, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(5, 0, 2, undefined, undefined, 0));
       faces.faces.push(new THREE.Face3(5, 2, 1, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(6, 4, 3, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(3, 4, 1, undefined, undefined, 0));
       faces.faces.push(new THREE.Face3(3, 1, 2, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(5, 6, 0, undefined, undefined, 0));
    
    
       faces.computeFaceNormals();
       faces.computeVertexNormals();
    
       var object = new THREE.Mesh(faces, faces_material);
       obj.add(object);
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 0, 1));
       allpoints.push(new PMPoint(1, 0, 0));
       allpoints.push(new PMPoint(0, 0, 0));
       allpoints.push(new PMPoint(0, 1, 0));
       allpoints.push(new PMPoint(1, 1, 0));
       allpoints.push(new PMPoint(1, 0, 1));
       allpoints.push(new PMPoint(0, 1, 1));
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       <!-- Edge style -->
       var line_material = new THREE.LineBasicMaterial ( {color: 0x000000, linewidth: 1.5, } );
    
       line_material.side = THREE.DoubleSide;
       line_material.transparent = true;
    
       <!-- EDGES -->
       var length = allpoints[5].vector.distanceTo(allpoints[0].vector)-allpoints[5].sphere.geometry.parameters.radius;
       var direction = allpoints[5].vector.clone()
       direction.sub(allpoints[0].vector)
       direction.normalize();
       var line = new THREE.ArrowHelper(direction, allpoints[0].vector, length, 0x000000, 0.2*length, 0.07);
    	obj.add(line);
       var length = allpoints[6].vector.distanceTo(allpoints[0].vector)-allpoints[6].sphere.geometry.parameters.radius;
       var direction = allpoints[6].vector.clone()
       direction.sub(allpoints[0].vector)
       direction.normalize();
       var line = new THREE.ArrowHelper(direction, allpoints[0].vector, length, 0x000000, 0.2*length, 0.07);
    	obj.add(line);
       var length = allpoints[4].vector.distanceTo(allpoints[1].vector)-allpoints[4].sphere.geometry.parameters.radius;
       var direction = allpoints[4].vector.clone()
       direction.sub(allpoints[1].vector)
       direction.normalize();
       var line = new THREE.ArrowHelper(direction, allpoints[1].vector, length, 0x000000, 0.2*length, 0.07);
    	obj.add(line);
       var length = allpoints[5].vector.distanceTo(allpoints[1].vector)-allpoints[5].sphere.geometry.parameters.radius;
       var direction = allpoints[5].vector.clone()
       direction.sub(allpoints[1].vector)
       direction.normalize();
       var line = new THREE.ArrowHelper(direction, allpoints[1].vector, length, 0x000000, 0.2*length, 0.07);
    	obj.add(line);
       var length = allpoints[0].vector.distanceTo(allpoints[2].vector)-allpoints[0].sphere.geometry.parameters.radius;
       var direction = allpoints[0].vector.clone()
       direction.sub(allpoints[2].vector)
       direction.normalize();
       var line = new THREE.ArrowHelper(direction, allpoints[2].vector, length, 0x000000, 0.2*length, 0.07);
    	obj.add(line);
       var length = allpoints[1].vector.distanceTo(allpoints[2].vector)-allpoints[1].sphere.geometry.parameters.radius;
       var direction = allpoints[1].vector.clone()
       direction.sub(allpoints[2].vector)
       direction.normalize();
       var line = new THREE.ArrowHelper(direction, allpoints[2].vector, length, 0x000000, 0.2*length, 0.07);
    	obj.add(line);
       var length = allpoints[3].vector.distanceTo(allpoints[2].vector)-allpoints[3].sphere.geometry.parameters.radius;
       var direction = allpoints[3].vector.clone()
       direction.sub(allpoints[2].vector)
       direction.normalize();
       var line = new THREE.ArrowHelper(direction, allpoints[2].vector, length, 0x000000, 0.2*length, 0.07);
    	obj.add(line);
       var length = allpoints[4].vector.distanceTo(allpoints[3].vector)-allpoints[4].sphere.geometry.parameters.radius;
       var direction = allpoints[4].vector.clone()
       direction.sub(allpoints[3].vector)
       direction.normalize();
       var line = new THREE.ArrowHelper(direction, allpoints[3].vector, length, 0x000000, 0.2*length, 0.07);
    	obj.add(line);
       var length = allpoints[6].vector.distanceTo(allpoints[3].vector)-allpoints[6].sphere.geometry.parameters.radius;
       var direction = allpoints[6].vector.clone()
       direction.sub(allpoints[3].vector)
       direction.normalize();
       var line = new THREE.ArrowHelper(direction, allpoints[3].vector, length, 0x000000, 0.2*length, 0.07);
    	obj.add(line);
       scene.add(obj);
       all_objects.push(obj);
    
    // COMMON_CODE_BLOCK_BEGIN
    var xRotationEnabled = false;
    var yRotationEnabled = false;
    var zRotationEnabled = false;
    var rotationSpeedFactor = 1;
    var settingsShown = false;
    var labelsShown = true;
    var intervals = [];
    var timeouts = [];
    var explodingSpeed = 0.05;
    var explodeScale = 0;
    var XMLS = new XMLSerializer();
    var svgElement;
    var renderId;
    
    	var render = function () {
    
    		renderId = requestAnimationFrame(render);
    
    //		comment in for automatic explosion
    //		explode(updateFactor());
    
    		var phi = 0.02 * rotationSpeedFactor;
    
    		if (xRotationEnabled){
    			scene.rotation.x += phi;
    		}
    		if(yRotationEnabled){
    			scene.rotation.y += phi;
    		}
    		if(zRotationEnabled){
    			scene.rotation.z += phi;
    		}
    
    		controls.update();
    		renderer.render(scene, camera);
    	};
    
    	render();
    
    	function computeCentroid(geom) {
    		centroid = new THREE.Vector3();
    		geom.vertices.forEach(function(v) {
    			centroid.add(v);			
    		});
    		centroid.divideScalar(geom.vertices.length);
    		return centroid;
    	}
    
    	function changeTransparency(event){
    		var opacity = 1-Number(event.currentTarget.value);
    		for (var i=0; i<all_objects.length; i++){
    			for (var j=0; j<all_objects[i].children.length; j++){
    				if (all_objects[i].children[j].material.type == "MultiMaterial") {
    					for (var k=0; k<all_objects[i].children[j].material.materials.length; k++){
    						all_objects[i].children[j].material.materials[k].opacity = opacity;
    						all_objects[i].children[j].material.materials[k].depthWrite = opacity < 0.5 ? false : true;
    						all_objects[i].children[j].material.materials[k].depthTest = opacity < 0.5 ? false : true;
    					}
    				} else if (all_objects[i].children[j].material.transparent && 
    							  all_objects[i].children[j].material.type == "MeshBasicMaterial" &&
    							  all_objects[i].children[j].geometry.type == "Geometry"){
    					all_objects[i].children[j].material.opacity = opacity;
    					all_objects[i].children[j].material.depthWrite = opacity < 0.5 ? false : true;
    					all_objects[i].children[j].material.depthTest = opacity < 0.5 ? false : true;
    				}
    			}
    		}
    	}
    
    	function changeRotationX(event){
    		xRotationEnabled = event.currentTarget.checked;
    	}	
    
    	function changeRotationY(event){
    		yRotationEnabled = event.currentTarget.checked;
    	}	
    
    	function changeRotationZ(event){
    		zRotationEnabled = event.currentTarget.checked;
    	}	
    
    
    	function changeRotationSpeedFactor(event){
    		rotationSpeedFactor = Number(event.currentTarget.value);
    	}
    
    	function resetScene(){
    		scene.rotation.set(0,0,0);
    		camera.position.set(0,0,5);
    		camera.up.set(0,1,0);
    	}
    
    	function showSettings(event){
    		event.currentTarget.style.display = 'none';
    		document.getElementById('settings_6').style.position = 'absolute';
    		document.getElementById('settings_6').style.display = 'block';
    		document.getElementById('showSettingsButton_6').style.display = 'none';
    		document.getElementById('hideSettingsButton_6').style.display = 'block';
    		settingsShown = true;
    	}
    
    	function hideSettings(event){
    		event.currentTarget.style.display = 'none';
    		document.getElementById('settings_6').style.display = 'none';
    		document.getElementById('hideSettingsButton_6').style.display = 'none';
    		document.getElementById('showSettingsButton_6').style.display = 'block';
    		settingsShown = false;
    	}
    
    
    
    	var pos = 150* Math.PI;
    
    	function updateFactor() {
    		pos++;
    		return Math.sin(.01*pos)+1;
    	}
    
    	function makelabel(message, x, y, z, params) {
    		var spritey = textSprite( message, params );
    		spritey.position.set(x, y, z);
    		obj.add(spritey);
    	}
    
    	function textSprite(message, parameters)
    	{
    		if ( parameters === undefined ) parameters = {};
    
    		var fontface = "Helvetica";
    
    		var fontsize = parameters.hasOwnProperty("fontsize") ? 
    			parameters["fontsize"] : 18;
    		fontsize = fontsize*10;
    
    		var canvas = document.createElement('canvas');
    		var size = 1024;
    		canvas.width = size;
    		canvas.height = size;
    		var context = canvas.getContext('2d');
    		context.font = fontsize + "px " + fontface;
    
    		// text color
    		context.fillStyle = "rgba(0, 0, 0, 1.0)";
    
    		context.fillText(message, size/2, size/2);
    
    		// canvas contents will be used for a texture
    		var texture = new THREE.Texture(canvas);
    		texture.needsUpdate = true;
    
    		var spriteMaterial = new THREE.SpriteMaterial(
    			{map: texture, useScreenCoordinates: false});
    		var sprite = new THREE.Sprite(spriteMaterial);
    		return sprite;
    	}
    
    	function takeSvgScreenshot(){
    		if (labelsShown){
    			hideLabels();
    		}
    		svgRenderer.render(scene,camera);
    		svgElement = XMLS.serializeToString(svgRenderer.domElement);
    		
    		if (labelsShown){
    			displayLabels();
    		}
    
    		if (document.getElementById('tab_6').checked){
    			//show in new tab
    			var myWindow = window.open("","");
    			myWindow.document.body.innerHTML = svgElement;
    		} else{
    			// download svg file 
    			download("screenshot.svg", svgElement);
    		}
    	}
    		
    
    	function showOrHideObject(event){
    		var nr = Number(event.currentTarget.name);
    		all_objects[nr].visible = event.currentTarget.checked;
    	}
    
    	function displayOrHideOptionsRecursive( obj ) {
    		for (var j=0; j<obj.children.length; j++) {
    			var child = obj.children[j];
    			if (child.material===undefined && child) {
    				displayOrHideOptionsRecursive( child );
    			} else {
    				if (child.material.type == "MultiMaterial") {
    					for (var k=0; k<child.material.materials.length; k++) {
    						if (child.material.materials[k].transparent) {
    							document.getElementById('transparency_6').style.display = 'block';
    							document.getElementById('transparencyRange_6').value = 1 - 
    								child.material.materials[k].opacity;
    							return;
    						}
    					}
    				} else if (	child.material.transparent && 
    								child.material.type == "MeshBasicMaterial" &&
    								child.geometry.type == "Geometry"){
    					document.getElementById('transparency_6').style.display = 'block';
    					return;
    				}
    			}
    		}
    	}
    
    	function displayOrHideOptions() {
    		for (var i=0; i<all_objects.length; i++) {
    			var obj = all_objects[i];
    			displayOrHideOptionsRecursive( obj );
    		}
    	}
    
    	displayOrHideOptions()
    
    
    
    
    // ---------------------- EXPLOSION ------------------------------------------------
    // ---------------------------------------------------------------------------------
    
    	function explode(factor) {
    		var obj, c;
    		var c0 = centroids[0];
    		for (var i = 0; i<centroids.length; ++i) {
    			c = centroids[i];
    			obj = all_objects[all_objects.length - centroids.length + i];
    			obj.position.set(c.x*factor, c.y*factor, c.z*factor);
    		}	
    	}
    
    	function triggerExplode(event){
    		explodeScale = Number(event.currentTarget.value);
    		explode(explodeScale);
    	}
    
    	function setExplodingSpeed(event){
    		explodingSpeed = Number(event.currentTarget.value);
    	}
    
    	function triggerAutomaticExplode(event){
    		if (event.currentTarget.checked){
    			startExploding();
    		} else {
    			clearIntervals();
    		}	
    	}
    
    	function startExploding(){
    		intervals.push(setInterval(explodingInterval, 25));
    	}
    
    
    	function explodingInterval(){
    		explodeScale += explodingSpeed;
    		if (explodeScale <= 6){ 
    			explode(explodeScale);
    		}
    		else{
    			explode(6);
    			explodeScale = 6;
    			clearIntervals();
    			timeouts.push(setTimeout(startUnexploding, 3000));
    		}
    		document.getElementById('explodeRange_6').value = explodeScale;
    	}
    
    
    	function startUnexploding(){
    		intervals.push(setInterval(unexplodingInterval, 25));
    	}
    
    	function unexplodingInterval(){
    		explodeScale -= explodingSpeed;
    		if (explodeScale >= 0){	
    			explode(explodeScale);
    		}
    		else {
    			explode(0);
    			explodeScale = 0;
    			clearIntervals();
    			timeouts.push(setTimeout(startExploding, 3000));
    		}
    		document.getElementById('explodeRange_6').value = explodeScale;
    	}
    
    	function clearIntervals(){
    		intervals.forEach(function(interval){
    			clearInterval(interval);
    		});
    		intervals = [];
    		timeouts.forEach(function(timeout){
    			clearTimeout(timeout);
    		});
    		timeouts = [];
    	}
    
    			
    
    	// append checkboxes for displaying or hiding objects
    	var shownObjectsList = document.getElementById('shownObjectsList_6');
    	for (var i=0; i<all_objects.length; i++){
    		var objNode = document.createElement('span');
    		objNode.innerHTML = objectnames[i] + '<br>';
    		var checkbox = document.createElement('input');
    		checkbox.type = 'checkbox';
    		checkbox.checked = true;
    		checkbox.name = String(i);
    		checkbox.onchange = showOrHideObject;
    		shownObjectsList.appendChild(checkbox);
    		shownObjectsList.appendChild(objNode);
    	}
    
    	function displayLabels(){
    		for (var i=0; i<all_objects.length; i++){
    			for (var j=0; j<all_objects[i].children.length; j++){
    				var child = all_objects[i].children[j];
    				if (child.type == 'Sprite'){
    					child.visible = true;
    				}
    			}
    		}
    	}
    
    	function hideLabels(){
    		for (var i=0; i<all_objects.length; i++){
    			for (var j=0; j<all_objects[i].children.length; j++){
    				var child = all_objects[i].children[j];
    				if (child.type == 'Sprite'){
    					child.visible = false;
    				}
    			}
    		}
    	}
    
    	function displayOrHideLabels(event){
    		if (event.currentTarget.checked){
    			displayLabels();
    			labelsShown = true;
    		} else {
    			hideLabels();
    			labelsShown = false;
    		}
    	}
    
    	function download(filename, text) {
    	  var element = document.createElement('a');
    	  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
    	  element.setAttribute('download', filename);
    
    	  element.style.display = 'none';
    	  document.body.appendChild(element);
    
    	  element.click();
    
    	  document.body.removeChild(element);
    	}
    
    var tempobj;
    tempobj = document.getElementById('explodeRange_6');
    if (tempobj) {
       tempobj.oninput = triggerExplode;
       document.getElementById('explodeCheckbox_6').onchange = triggerAutomaticExplode;
       document.getElementById('explodingSpeedRange_6').oninput = setExplodingSpeed;
    }
    tempobj = document.getElementById('foldRange_6');
    if (tempobj) {
       tempobj.oninput = fold;
    }
    document.getElementById('transparencyRange_6').oninput = changeTransparency;
    document.getElementById('changeRotationX_6').onchange = changeRotationX;
    document.getElementById('changeRotationY_6').onchange = changeRotationY;
    document.getElementById('changeRotationZ_6').onchange = changeRotationZ;
    document.getElementById('resetButton_6').onclick = resetScene;
    document.getElementById('rotationSpeedRange_6').oninput = changeRotationSpeedFactor;
    document.getElementById('labelsCheckboxInput_6').onchange = displayOrHideLabels;
    document.getElementById('takeScreenshot_6').onclick = takeSvgScreenshot;
    document.getElementById('showSettingsButton_6').onclick = showSettings;
    document.getElementById('hideSettingsButton_6').onclick = hideSettings;
    
    	
    
    // ------------------ SHORTCUTS --------------------------------------------
    // -------------------------------------------------------------------------
    
    /**
     * http://www.openjs.com/scripts/events/keyboard_shortcuts/
     * Version : 2.01.B
     * By Binny V A
     * License : BSD
     */
    shortcut = {
    	'all_shortcuts':{},//All the shortcuts are stored in this array
    	'add': function(shortcut_combination,callback,opt) {
    		//Provide a set of default options
    		var default_options = {
    			'type':'keydown',
    			'propagate':false,
    			'disable_in_input':false,
    			'target':document,
    			'keycode':false
    		}
    		if(!opt) opt = default_options;
    		else {
    			for(var dfo in default_options) {
    				if(typeof opt[dfo] == 'undefined') opt[dfo] = default_options[dfo];
    			}
    		}
    
    		var ele = opt.target;
    		if(typeof opt.target == 'string') ele = document.getElementById(opt.target);
    		var ths = this;
    		shortcut_combination = shortcut_combination.toLowerCase();
    
    		//The function to be called at keypress
    		var func = function(e) {
    			e = e || window.event;
    			
    			if(opt['disable_in_input']) { //Don't enable shortcut keys in Input, Textarea fields
    				var element;
    				if(e.target) element=e.target;
    				else if(e.srcElement) element=e.srcElement;
    				if(element.nodeType==3) element=element.parentNode;
    
    				if(element.tagName == 'INPUT' || element.tagName == 'TEXTAREA') return;
    			}
    	
    			//Find Which key is pressed
    			if (e.keyCode) code = e.keyCode;
    			else if (e.which) code = e.which;
    			var character = String.fromCharCode(code).toLowerCase();
    			
    			if(code == 188) character=","; //If the user presses , when the type is onkeydown
    			if(code == 190) character="."; //If the user presses , when the type is onkeydown
    
    			var keys = shortcut_combination.split("+");
    			//Key Pressed - counts the number of valid keypresses - if it is same as the number of keys, the shortcut function is invoked
    			var kp = 0;
    			
    			//Work around for stupid Shift key bug created by using lowercase - as a result the shift+num combination was broken
    			var shift_nums = {
    				"`":"~",
    				"1":"!",
    				"2":"@",
    				"3":"#",
    				"4":"$",
    				"5":"%",
    				"6":"^",
    				"7":"&",
    				"8":"*",
    				"9":"(",
    				"0":")",
    				"-":"_",
    				"=":"+",
    				";":":",
    				"'":"\"",
    				",":"<",
    				".":">",
    				"/":"?",
    				"\\":"|"
    			}
    			//Special Keys - and their codes
    			var special_keys = {
    				'esc':27,
    				'escape':27,
    				'tab':9,
    				'space':32,
    				'return':13,
    				'enter':13,
    				'backspace':8,
    	
    				'scrolllock':145,
    				'scroll_lock':145,
    				'scroll':145,
    				'capslock':20,
    				'caps_lock':20,
    				'caps':20,
    				'numlock':144,
    				'num_lock':144,
    				'num':144,
    				
    				'pause':19,
    				'break':19,
    				
    				'insert':45,
    				'home':36,
    				'delete':46,
    				'end':35,
    				
    				'pageup':33,
    				'page_up':33,
    				'pu':33,
    	
    				'pagedown':34,
    				'page_down':34,
    				'pd':34,
    	
    				'left':37,
    				'up':38,
    				'right':39,
    				'down':40,
    	
    				'f1':112,
    				'f2':113,
    				'f3':114,
    				'f4':115,
    				'f5':116,
    				'f6':117,
    				'f7':118,
    				'f8':119,
    				'f9':120,
    				'f10':121,
    				'f11':122,
    				'f12':123
    			}
    	
    			var modifiers = { 
    				shift: { wanted:false, pressed:false},
    				ctrl : { wanted:false, pressed:false},
    				alt  : { wanted:false, pressed:false},
    				meta : { wanted:false, pressed:false}	//Meta is Mac specific
    			};
                            
    			if(e.ctrlKey)	modifiers.ctrl.pressed = true;
    			if(e.shiftKey)	modifiers.shift.pressed = true;
    			if(e.altKey)	modifiers.alt.pressed = true;
    			if(e.metaKey)   modifiers.meta.pressed = true;
                            
    			for(var i=0; k=keys[i],i<keys.length; i++) {
    				//Modifiers
    				if(k == 'ctrl' || k == 'control') {
    					kp++;
    					modifiers.ctrl.wanted = true;
    
    				} else if(k == 'shift') {
    					kp++;
    					modifiers.shift.wanted = true;
    
    				} else if(k == 'alt') {
    					kp++;
    					modifiers.alt.wanted = true;
    				} else if(k == 'meta') {
    					kp++;
    					modifiers.meta.wanted = true;
    				} else if(k.length > 1) { //If it is a special key
    					if(special_keys[k] == code) kp++;
    					
    				} else if(opt['keycode']) {
    					if(opt['keycode'] == code) kp++;
    
    				} else { //The special keys did not match
    					if(character == k) kp++;
    					else {
    						if(shift_nums[character] && e.shiftKey) { //Stupid Shift key bug created by using lowercase
    							character = shift_nums[character]; 
    							if(character == k) kp++;
    						}
    					}
    				}
    			}
    			
    			if(kp == keys.length && 
    						modifiers.ctrl.pressed == modifiers.ctrl.wanted &&
    						modifiers.shift.pressed == modifiers.shift.wanted &&
    						modifiers.alt.pressed == modifiers.alt.wanted &&
    						modifiers.meta.pressed == modifiers.meta.wanted) {
    				callback(e);
    	
    				if(!opt['propagate']) { //Stop the event
    					//e.cancelBubble is supported by IE - this will kill the bubbling process.
    					e.cancelBubble = true;
    					e.returnValue = false;
    	
    					//e.stopPropagation works in Firefox.
    					if (e.stopPropagation) {
    						e.stopPropagation();
    						e.preventDefault();
    					}
    					return false;
    				}
    			}
    		}
    		this.all_shortcuts[shortcut_combination] = {
    			'callback':func, 
    			'target':ele, 
    			'event': opt['type']
    		};
    		//Attach the function with the event
    		if(ele.addEventListener) ele.addEventListener(opt['type'], func, false);
    		else if(ele.attachEvent) ele.attachEvent('on'+opt['type'], func);
    		else ele['on'+opt['type']] = func;
    	},
    
    	//Remove the shortcut - just specify the shortcut and I will remove the binding
    	'remove':function(shortcut_combination) {
    		shortcut_combination = shortcut_combination.toLowerCase();
    		var binding = this.all_shortcuts[shortcut_combination];
    		delete(this.all_shortcuts[shortcut_combination])
    		if(!binding) return;
    		var type = binding['event'];
    		var ele = binding['target'];
    		var callback = binding['callback'];
    
    		if(ele.detachEvent) ele.detachEvent('on'+type, callback);
    		else if(ele.removeEventListener) ele.removeEventListener(type, callback, false);
    		else ele['on'+type] = false;
    	}
    }
    
    shortcut.add("Alt+Left",function() {
    	var event = new Event('click');
    	if (settingsShown){
    		document.getElementById('hideSettingsButton_6').dispatchEvent(event);
    	} else{
    		document.getElementById('showSettingsButton_6').dispatchEvent(event);
    	}
    });
    
    if (foldable) moveToBaryCenter();
    
    
    });});
    // COMMON_CODE_BLOCK_END
    </script>
    
    </body>
    </html>



.. figure:: attachment:c3t_graph.gif
   :alt: {{ :tutorial:c3t_graph.gif?300 }}

   {{ :tutorial:c3t_graph.gif?300 }}

The minimal and maximal faces can be visualized via


::

    polymake> $p->VISUAL->MIN_MAX_FACE;


.. raw:: html

    <!--
    polymake for andrew
    Thu Mar 28 11:14:58 2019
    c3t
    -->
    
    
    <html>
       <head>
          <title>c3t</title>
          <style>
    /*
    // COMMON_CODE_BLOCK_BEGIN
    */
             html{overflow: scroll;}
             body { font-family: Arial, Helvetica, sans-serif}
             strong{font-size: 18px;}
             canvas { z-index: 8; }
             input[type='range'] {}
             input[type='radio'] {margin-left:0;}
             input[type='checkbox'] {margin-right:7px; margin-left: 0px; padding-left:0px;}
             .group{padding-bottom: 40px;}
             .settings * {z-index: 11; }
             .settings{z-index: 10; margin-left: 30px; display: none; width: 14em; height: 90%; border: solid 1px silver; padding: 2px; overflow-y: scroll; background-color: white }
             .indented{margin-left: 20px; margin-top: 15px; padding-bottom: 0px;} 
             .shownObjectsList{overflow: auto; max-width: 150px; max-height: 150px;}
             .showSettingsButton{display: block; z-index: 12; position: absolute }
             .hideSettingsButton{display: none; z-index: 12; position: absolute; opacity: 0.5}
             .resetButton{margin-top: 20px;}
             button{margin-left: 0;}
             img{cursor: pointer;}
             .suboption{padding-top: 30px;}
             .transparency{display: none;}
             .labelsCheckbox{margin-top: 10px;}
    
    
             input[type=range] {
               -webkit-appearance: none;
               padding:0; 
               width:90%; 
               margin-left: auto;
               margin-right: auto;
               margin-top: 20px;
               display: block;	
             }
             input[type=range]:focus {
               outline: none;
             }
             input[type=range]::-webkit-slider-runnable-track {
               height: 4px;
               cursor: pointer;
               animate: 0.2s;
               box-shadow: 0px 0px 0px #000000;
               background: #E3E3E3;
               border-radius: 0px;
               border: 0px solid #000000;
             }
             input[type=range]::-webkit-slider-thumb {
               box-shadow: 1px 1px 2px #B8B8B8;
               border: 1px solid #ABABAB;
               height: 13px;
               width: 25px;
               border-radius: 20px;
               background: #E0E0E0;
               cursor: pointer;
               -webkit-appearance: none;
               margin-top: -5px;
             }
             input[type=range]:focus::-webkit-slider-runnable-track {
               background: #E3E3E3;
             }
             input[type=range]::-moz-range-track {
               height: 4px;
               cursor: pointer;
               animate: 0.2s;
               box-shadow: 0px 0px 0px #000000;
               background: #E3E3E3;
               border-radius: 0px;
               border: 0px solid #000000;
             }
             input[type=range]::-moz-range-thumb {
               box-shadow: 1px 1px 2px #B8B8B8;
               border: 1px solid #ABABAB;
               height: 13px;
               width: 25px;
               border-radius: 20px;
               background: #E0E0E0;
               cursor: pointer;
             }
             input[type=range]::-ms-track {
               height: 4px;
               cursor: pointer;
               animate: 0.2s;
               background: transparent;
               border-color: transparent;
               color: transparent;
             }
             input[type=range]::-ms-fill-lower {
               background: #E3E3E3;
               border: 0px solid #000000;
               border-radius: 0px;
               box-shadow: 0px 0px 0px #000000;
             }
             input[type=range]::-ms-fill-upper {
               background: #E3E3E3;
               border: 0px solid #000000;
               border-radius: 0px;
               box-shadow: 0px 0px 0px #000000;
             }
             input[type=range]::-ms-thumb {
               box-shadow: 1px 1px 2px #B8B8B8;
               border: 1px solid #ABABAB;
               height: 13px;
               width: 25px;
               border-radius: 20px;
               background: #E0E0E0;
               cursor: pointer;
             }
             input[type=range]:focus::-ms-fill-lower {
               background: #E3E3E3;
             }
             input[type=range]:focus::-ms-fill-upper {
               background: #E3E3E3;
             }
    /*
    // COMMON_CODE_BLOCK_END
    */
    		</style>
       </head>
    
    <body>
    
    		<div id='settings_7' class='settings'>
    			<div class=group id='transparency_7' class='transparency'>
    				<strong>Transparency</strong>
    				<input id='transparencyRange_7' type='range' min=0 max=1 step=0.01 value=0>
    			</div>
    			
    			<div class=group id='rotation_7'>
    				<strong>Rotation</strong>
    				<div class=indented>
    					<div><input type='checkbox' id='changeRotationX_7'> x-axis</div>
    					<div><input type='checkbox' id='changeRotationY_7'> y-axis</div>
    					<div><input type='checkbox' id='changeRotationZ_7'> z-axis</div>
    					<button id='resetButton_7' class='resetButton' >Reset</button>
    				</div>
    
    				<div class=suboption>Rotation speed</div>
    				<input id='rotationSpeedRange_7' type='range' min=0 max=5 step=0.01 value=2>
    
    			</div>
    
    
    			<div class=group id='display_7'>
    				<strong>Display</strong>
    				<div class=indented>
    					<div id='shownObjectsList_7' class='shownObjectsList'></div>
    					<div class='labelsCheckbox'><input type='checkbox' id='labelsCheckboxInput_7' checked>Labels</div>
    				</div>
    			</div>
    
    
    			<div class=group id='svg_7'>
    				<strong>SVG</strong>
    				<div class=indented>
    					<form>
    						<input type="radio" name='screenshotMode' value='download' id='download_7' checked> Download<br>
    						<input type="radio" name='screenshotMode' value='tab' id='tab_7' > New tab<br>
    					</form>
    					<button id='takeScreenshot_7'>Screenshot</button>
    				</div>
    			</div>
    
    		</div>	<!-- end of settings -->
    		<img id='hideSettingsButton_7' style="display: none" class='hideSettingsButton' src='/kernelspecs/polymake/close.svg' width=20px">
    		<img id='showSettingsButton_7' class='showSettingsButton' src='/kernelspecs/polymake/menu.svg' width=20px">
    <div id="model6193174868"></div>
    
    <script>
    requirejs.config({
      paths: {
        three: '/kernelspecs/polymake/three',
        Detector: '/kernelspecs/polymake/Detector',
        SVGRenderer: '/kernelspecs/polymake/SVGRenderer',
        CanvasRenderer: '/kernelspecs/polymake/CanvasRenderer',
        Projector: '/kernelspecs/polymake/Projector',
        TrackballControls: '/kernelspecs/polymake/TrackballControls'
      },
      shim: {
        'three':
        {
          exports: 'THREE'
        },
        'Detector':
        {
          deps: [ 'three' ],
          exports: 'Detector'
        },
        'SVGRenderer':
        {
          deps: [ 'three' ],
          exports: 'THREE.SVGRenderer'
        },
        'CanvasRenderer':
        {
          deps: [ 'three' ],
          exports: 'THREE.CanvasRenderer'
        },
        'Projector':
        {
          deps: [ 'three' ],
          exports: 'THREE.Projector'
        },
        'TrackballControls':
        {
          deps: [ 'three' ],
          exports: 'THREE.TrackballControls'
        }
      }
    });
    require(['three'],function(THREE){
        window.THREE = THREE;
      require(['Detector','SVGRenderer','CanvasRenderer','Projector','TrackballControls'],function(Detector,SVGRenderer,CanvasRenderer,Projector,TrackballControls){
          THREE.SVGRenderer = SVGRenderer;
          THREE.CanvasRenderer = CanvasRenderer;
          THREE.Projector = Projector;
          THREE.TrackballControls = TrackballControls;
    
    // COMMON_CODE_BLOCK_BEGIN
    	var foldable = false;
       var container = document.getElementById( 'model6193174868' );
       var renderer = Detector.webgl? new THREE.WebGLRenderer({antialias: true}): new THREE.CanvasRenderer({antialias: true});
    	var svgRenderer = new THREE.SVGRenderer({antialias: true});
                var box = document.getElementsByClassName( 'output_subarea' )[0];
             var notebook = document.getElementById( 'notebook_panel' );
    
       var width = box.clientWidth - 25;
       var height = notebook.clientHeight * 0.8;
       renderer.setSize(width, height);
       svgRenderer.setSize(width, height);
       renderer.setClearColor(0xFFFFFF, 1);
       svgRenderer.setClearColor(0xFFFFFF, 1);
    
       container.appendChild(renderer.domElement);
    
       var scene = new THREE.Scene();
       var camera = new THREE.PerspectiveCamera(75, width/height, 0.1, 1000);
    
       var renderid;
    
       camera.position.set(0, 0, 5);
       camera.lookAt(0, 0, 0);
       camera.up.set(0, 1, 0);
    
       // class to allow move points together with labels and spheres
       var PMPoint = function (x,y,z) {
          this.vector = new THREE.Vector3(x,y,z);
          this.sprite = null;
          this.sphere = null;
       }
       PMPoint.prototype.makelabel = function(label) {
          this.sprite = textSprite( label );
          this.sprite.position.copy(this.vector);
       }
       PMPoint.prototype.makesphere = function(radius,material) {
          this.sphere = new THREE.Mesh(new THREE.SphereGeometry(radius), material);
          this.sphere.position.copy(this.vector);
       }
    
       PMPoint.prototype.setX = function(x) {
          this.vector.setX(x);
          if (this.sprite) {
             this.sprite.position.setX(x);
          }
          if (this.sphere) {
             this.sphere.position.setX(x);
          }
       };
       PMPoint.prototype.setY = function(y) {
          this.vector.setY(y);
          if (this.sprite) {
             this.sprite.position.setY(y);
          }
          if (this.sphere) {
             this.sphere.position.setY(y);
          }
       };
       PMPoint.prototype.setZ = function(z) {
          this.vector.setZ(z);
          if (this.sprite) {
             this.sprite.position.setZ(z);
          }
          if (this.sphere) {
             this.sphere.position.setZ(z);
          }
       };
       PMPoint.prototype.set = function(x,y,z) {
          this.vector.set(x,y,z);
          if (this.sprite) {
             this.sprite.position.set(x,y,z);
          }
          if (this.sphere) {
             this.sphere.position.set(x,y,z);
          }
       };
       PMPoint.prototype.add = function(o) {
          if (this.sprite) {
             o.add(this.sprite);
          }
          if (this.sphere) {
             o.add(this.sphere);
          }
       };
    
    
       var controls = new THREE.TrackballControls(camera, container);
    	controls.zoomSpeed = 0.2;
    	controls.rotateSpeed = 4;
    
       var all_objects = [];
       var centroids = [];
       // select the target node
       var target = document.querySelector('#model6193174868');
    
       // create an observer instance
       var observer = new MutationObserver(function(mutations) {
          mutations.forEach(function(mutation) {
             if (mutation.removedNodes && mutation.removedNodes.length > 0) {
                cancelAnimationFrame(renderId);
                observer.disconnect();
                console.log("cancelled frame "+renderId);
             }
          });
       });
    
       // configuration of the observer:
       var config = { childList: true, characterData: true }
    
       // pass in the target node, as well as the observer options
       while (target) {
          if (target.className=="output") {
             observer.observe(target, config);
             break;
          }
          target = target.parentNode;
       }
    
    // COMMON_CODE_BLOCK_END
    
       var objectnames = ["c3t"];
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 0, 1));
       allpoints.push(new PMPoint(1, 0, 0));
       allpoints.push(new PMPoint(0, 0, 0));
       allpoints.push(new PMPoint(0, 1, 0));
       allpoints.push(new PMPoint(1, 1, 0));
       allpoints.push(new PMPoint(1, 0, 1));
       allpoints.push(new PMPoint(0, 1, 1));
    
       <!-- Vertex style -->
       var materials = [
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFFFF00, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
       ];
       for (index = 0; index < materials.length; ++index) {
          materials[index].side = THREE.DoubleSide;
       }
       var points_material = new THREE.MeshFaceMaterial ( materials );
    
    
       <!-- POINTS -->
       allpoints[0].makesphere(0.02,materials[0]);
       allpoints[1].makesphere(0.02,materials[1]);
       allpoints[2].makesphere(0.02,materials[2]);
       allpoints[3].makesphere(0.02,materials[3]);
       allpoints[4].makesphere(0.02,materials[4]);
       allpoints[5].makesphere(0.02,materials[5]);
       allpoints[6].makesphere(0.02,materials[6]);
       allpoints[0].makelabel("0");
       allpoints[1].makelabel("1");
       allpoints[2].makelabel("2");
       allpoints[3].makelabel("3");
       allpoints[4].makelabel("4");
       allpoints[5].makelabel("5");
       allpoints[6].makelabel("6");
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
       faces.vertices.push(allpoints[3].vector);
       faces.vertices.push(allpoints[4].vector);
       faces.vertices.push(allpoints[5].vector);
       faces.vertices.push(allpoints[6].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var materials = [
          new THREE.MeshBasicMaterial({ transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, color: 0x77EC9E, }),
          new THREE.MeshBasicMaterial({ transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, color: 0x77EC9E, }),
          new THREE.MeshBasicMaterial({ transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, color: 0x77EC9E, }),
          new THREE.MeshBasicMaterial({ transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, color: 0x77EC9E, }),
          new THREE.MeshBasicMaterial({ transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, color: 0x77EC9E, }),
          new THREE.MeshBasicMaterial({ transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, color: 0x77EC9E, }),
       ];
       for (index = 0; index < materials.length; ++index) {
          materials[index].side = THREE.DoubleSide;
       }
       var faces_material = new THREE.MeshFaceMaterial ( materials );
    
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(5, 4, 6, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(6, 3, 2, undefined, undefined, 1));
       faces.faces.push(new THREE.Face3(6, 2, 0, undefined, undefined, 1));
    
       faces.faces.push(new THREE.Face3(1, 4, 5, undefined, undefined, 2));
    
       faces.faces.push(new THREE.Face3(5, 0, 2, undefined, undefined, 3));
       faces.faces.push(new THREE.Face3(5, 2, 1, undefined, undefined, 3));
    
       faces.faces.push(new THREE.Face3(6, 4, 3, undefined, undefined, 4));
    
       faces.faces.push(new THREE.Face3(3, 4, 1, undefined, undefined, 5));
       faces.faces.push(new THREE.Face3(3, 1, 2, undefined, undefined, 5));
    
       faces.faces.push(new THREE.Face3(5, 6, 0, undefined, undefined, 6));
    
    
       faces.computeFaceNormals();
       faces.computeVertexNormals();
    
       var object = new THREE.Mesh(faces, faces_material);
       obj.add(object);
    
       <!-- Edge style -->
       var line_material = new THREE.LineBasicMaterial ( {color: 0x000000, linewidth: 1.5, } );
    
       line_material.side = THREE.DoubleSide;
       line_material.transparent = true;
    
       <!-- EDGES --> 
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[5].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[6].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[5].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[6].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[3].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[5].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
    // COMMON_CODE_BLOCK_BEGIN
    var xRotationEnabled = false;
    var yRotationEnabled = false;
    var zRotationEnabled = false;
    var rotationSpeedFactor = 1;
    var settingsShown = false;
    var labelsShown = true;
    var intervals = [];
    var timeouts = [];
    var explodingSpeed = 0.05;
    var explodeScale = 0;
    var XMLS = new XMLSerializer();
    var svgElement;
    var renderId;
    
    	var render = function () {
    
    		renderId = requestAnimationFrame(render);
    
    //		comment in for automatic explosion
    //		explode(updateFactor());
    
    		var phi = 0.02 * rotationSpeedFactor;
    
    		if (xRotationEnabled){
    			scene.rotation.x += phi;
    		}
    		if(yRotationEnabled){
    			scene.rotation.y += phi;
    		}
    		if(zRotationEnabled){
    			scene.rotation.z += phi;
    		}
    
    		controls.update();
    		renderer.render(scene, camera);
    	};
    
    	render();
    
    	function computeCentroid(geom) {
    		centroid = new THREE.Vector3();
    		geom.vertices.forEach(function(v) {
    			centroid.add(v);			
    		});
    		centroid.divideScalar(geom.vertices.length);
    		return centroid;
    	}
    
    	function changeTransparency(event){
    		var opacity = 1-Number(event.currentTarget.value);
    		for (var i=0; i<all_objects.length; i++){
    			for (var j=0; j<all_objects[i].children.length; j++){
    				if (all_objects[i].children[j].material.type == "MultiMaterial") {
    					for (var k=0; k<all_objects[i].children[j].material.materials.length; k++){
    						all_objects[i].children[j].material.materials[k].opacity = opacity;
    						all_objects[i].children[j].material.materials[k].depthWrite = opacity < 0.5 ? false : true;
    						all_objects[i].children[j].material.materials[k].depthTest = opacity < 0.5 ? false : true;
    					}
    				} else if (all_objects[i].children[j].material.transparent && 
    							  all_objects[i].children[j].material.type == "MeshBasicMaterial" &&
    							  all_objects[i].children[j].geometry.type == "Geometry"){
    					all_objects[i].children[j].material.opacity = opacity;
    					all_objects[i].children[j].material.depthWrite = opacity < 0.5 ? false : true;
    					all_objects[i].children[j].material.depthTest = opacity < 0.5 ? false : true;
    				}
    			}
    		}
    	}
    
    	function changeRotationX(event){
    		xRotationEnabled = event.currentTarget.checked;
    	}	
    
    	function changeRotationY(event){
    		yRotationEnabled = event.currentTarget.checked;
    	}	
    
    	function changeRotationZ(event){
    		zRotationEnabled = event.currentTarget.checked;
    	}	
    
    
    	function changeRotationSpeedFactor(event){
    		rotationSpeedFactor = Number(event.currentTarget.value);
    	}
    
    	function resetScene(){
    		scene.rotation.set(0,0,0);
    		camera.position.set(0,0,5);
    		camera.up.set(0,1,0);
    	}
    
    	function showSettings(event){
    		event.currentTarget.style.display = 'none';
    		document.getElementById('settings_7').style.position = 'absolute';
    		document.getElementById('settings_7').style.display = 'block';
    		document.getElementById('showSettingsButton_7').style.display = 'none';
    		document.getElementById('hideSettingsButton_7').style.display = 'block';
    		settingsShown = true;
    	}
    
    	function hideSettings(event){
    		event.currentTarget.style.display = 'none';
    		document.getElementById('settings_7').style.display = 'none';
    		document.getElementById('hideSettingsButton_7').style.display = 'none';
    		document.getElementById('showSettingsButton_7').style.display = 'block';
    		settingsShown = false;
    	}
    
    
    
    	var pos = 150* Math.PI;
    
    	function updateFactor() {
    		pos++;
    		return Math.sin(.01*pos)+1;
    	}
    
    	function makelabel(message, x, y, z, params) {
    		var spritey = textSprite( message, params );
    		spritey.position.set(x, y, z);
    		obj.add(spritey);
    	}
    
    	function textSprite(message, parameters)
    	{
    		if ( parameters === undefined ) parameters = {};
    
    		var fontface = "Helvetica";
    
    		var fontsize = parameters.hasOwnProperty("fontsize") ? 
    			parameters["fontsize"] : 18;
    		fontsize = fontsize*10;
    
    		var canvas = document.createElement('canvas');
    		var size = 1024;
    		canvas.width = size;
    		canvas.height = size;
    		var context = canvas.getContext('2d');
    		context.font = fontsize + "px " + fontface;
    
    		// text color
    		context.fillStyle = "rgba(0, 0, 0, 1.0)";
    
    		context.fillText(message, size/2, size/2);
    
    		// canvas contents will be used for a texture
    		var texture = new THREE.Texture(canvas);
    		texture.needsUpdate = true;
    
    		var spriteMaterial = new THREE.SpriteMaterial(
    			{map: texture, useScreenCoordinates: false});
    		var sprite = new THREE.Sprite(spriteMaterial);
    		return sprite;
    	}
    
    	function takeSvgScreenshot(){
    		if (labelsShown){
    			hideLabels();
    		}
    		svgRenderer.render(scene,camera);
    		svgElement = XMLS.serializeToString(svgRenderer.domElement);
    		
    		if (labelsShown){
    			displayLabels();
    		}
    
    		if (document.getElementById('tab_7').checked){
    			//show in new tab
    			var myWindow = window.open("","");
    			myWindow.document.body.innerHTML = svgElement;
    		} else{
    			// download svg file 
    			download("screenshot.svg", svgElement);
    		}
    	}
    		
    
    	function showOrHideObject(event){
    		var nr = Number(event.currentTarget.name);
    		all_objects[nr].visible = event.currentTarget.checked;
    	}
    
    	function displayOrHideOptionsRecursive( obj ) {
    		for (var j=0; j<obj.children.length; j++) {
    			var child = obj.children[j];
    			if (child.material===undefined && child) {
    				displayOrHideOptionsRecursive( child );
    			} else {
    				if (child.material.type == "MultiMaterial") {
    					for (var k=0; k<child.material.materials.length; k++) {
    						if (child.material.materials[k].transparent) {
    							document.getElementById('transparency_7').style.display = 'block';
    							document.getElementById('transparencyRange_7').value = 1 - 
    								child.material.materials[k].opacity;
    							return;
    						}
    					}
    				} else if (	child.material.transparent && 
    								child.material.type == "MeshBasicMaterial" &&
    								child.geometry.type == "Geometry"){
    					document.getElementById('transparency_7').style.display = 'block';
    					return;
    				}
    			}
    		}
    	}
    
    	function displayOrHideOptions() {
    		for (var i=0; i<all_objects.length; i++) {
    			var obj = all_objects[i];
    			displayOrHideOptionsRecursive( obj );
    		}
    	}
    
    	displayOrHideOptions()
    
    
    
    
    // ---------------------- EXPLOSION ------------------------------------------------
    // ---------------------------------------------------------------------------------
    
    	function explode(factor) {
    		var obj, c;
    		var c0 = centroids[0];
    		for (var i = 0; i<centroids.length; ++i) {
    			c = centroids[i];
    			obj = all_objects[all_objects.length - centroids.length + i];
    			obj.position.set(c.x*factor, c.y*factor, c.z*factor);
    		}	
    	}
    
    	function triggerExplode(event){
    		explodeScale = Number(event.currentTarget.value);
    		explode(explodeScale);
    	}
    
    	function setExplodingSpeed(event){
    		explodingSpeed = Number(event.currentTarget.value);
    	}
    
    	function triggerAutomaticExplode(event){
    		if (event.currentTarget.checked){
    			startExploding();
    		} else {
    			clearIntervals();
    		}	
    	}
    
    	function startExploding(){
    		intervals.push(setInterval(explodingInterval, 25));
    	}
    
    
    	function explodingInterval(){
    		explodeScale += explodingSpeed;
    		if (explodeScale <= 6){ 
    			explode(explodeScale);
    		}
    		else{
    			explode(6);
    			explodeScale = 6;
    			clearIntervals();
    			timeouts.push(setTimeout(startUnexploding, 3000));
    		}
    		document.getElementById('explodeRange_7').value = explodeScale;
    	}
    
    
    	function startUnexploding(){
    		intervals.push(setInterval(unexplodingInterval, 25));
    	}
    
    	function unexplodingInterval(){
    		explodeScale -= explodingSpeed;
    		if (explodeScale >= 0){	
    			explode(explodeScale);
    		}
    		else {
    			explode(0);
    			explodeScale = 0;
    			clearIntervals();
    			timeouts.push(setTimeout(startExploding, 3000));
    		}
    		document.getElementById('explodeRange_7').value = explodeScale;
    	}
    
    	function clearIntervals(){
    		intervals.forEach(function(interval){
    			clearInterval(interval);
    		});
    		intervals = [];
    		timeouts.forEach(function(timeout){
    			clearTimeout(timeout);
    		});
    		timeouts = [];
    	}
    
    			
    
    	// append checkboxes for displaying or hiding objects
    	var shownObjectsList = document.getElementById('shownObjectsList_7');
    	for (var i=0; i<all_objects.length; i++){
    		var objNode = document.createElement('span');
    		objNode.innerHTML = objectnames[i] + '<br>';
    		var checkbox = document.createElement('input');
    		checkbox.type = 'checkbox';
    		checkbox.checked = true;
    		checkbox.name = String(i);
    		checkbox.onchange = showOrHideObject;
    		shownObjectsList.appendChild(checkbox);
    		shownObjectsList.appendChild(objNode);
    	}
    
    	function displayLabels(){
    		for (var i=0; i<all_objects.length; i++){
    			for (var j=0; j<all_objects[i].children.length; j++){
    				var child = all_objects[i].children[j];
    				if (child.type == 'Sprite'){
    					child.visible = true;
    				}
    			}
    		}
    	}
    
    	function hideLabels(){
    		for (var i=0; i<all_objects.length; i++){
    			for (var j=0; j<all_objects[i].children.length; j++){
    				var child = all_objects[i].children[j];
    				if (child.type == 'Sprite'){
    					child.visible = false;
    				}
    			}
    		}
    	}
    
    	function displayOrHideLabels(event){
    		if (event.currentTarget.checked){
    			displayLabels();
    			labelsShown = true;
    		} else {
    			hideLabels();
    			labelsShown = false;
    		}
    	}
    
    	function download(filename, text) {
    	  var element = document.createElement('a');
    	  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
    	  element.setAttribute('download', filename);
    
    	  element.style.display = 'none';
    	  document.body.appendChild(element);
    
    	  element.click();
    
    	  document.body.removeChild(element);
    	}
    
    var tempobj;
    tempobj = document.getElementById('explodeRange_7');
    if (tempobj) {
       tempobj.oninput = triggerExplode;
       document.getElementById('explodeCheckbox_7').onchange = triggerAutomaticExplode;
       document.getElementById('explodingSpeedRange_7').oninput = setExplodingSpeed;
    }
    tempobj = document.getElementById('foldRange_7');
    if (tempobj) {
       tempobj.oninput = fold;
    }
    document.getElementById('transparencyRange_7').oninput = changeTransparency;
    document.getElementById('changeRotationX_7').onchange = changeRotationX;
    document.getElementById('changeRotationY_7').onchange = changeRotationY;
    document.getElementById('changeRotationZ_7').onchange = changeRotationZ;
    document.getElementById('resetButton_7').onclick = resetScene;
    document.getElementById('rotationSpeedRange_7').oninput = changeRotationSpeedFactor;
    document.getElementById('labelsCheckboxInput_7').onchange = displayOrHideLabels;
    document.getElementById('takeScreenshot_7').onclick = takeSvgScreenshot;
    document.getElementById('showSettingsButton_7').onclick = showSettings;
    document.getElementById('hideSettingsButton_7').onclick = hideSettings;
    
    	
    
    // ------------------ SHORTCUTS --------------------------------------------
    // -------------------------------------------------------------------------
    
    /**
     * http://www.openjs.com/scripts/events/keyboard_shortcuts/
     * Version : 2.01.B
     * By Binny V A
     * License : BSD
     */
    shortcut = {
    	'all_shortcuts':{},//All the shortcuts are stored in this array
    	'add': function(shortcut_combination,callback,opt) {
    		//Provide a set of default options
    		var default_options = {
    			'type':'keydown',
    			'propagate':false,
    			'disable_in_input':false,
    			'target':document,
    			'keycode':false
    		}
    		if(!opt) opt = default_options;
    		else {
    			for(var dfo in default_options) {
    				if(typeof opt[dfo] == 'undefined') opt[dfo] = default_options[dfo];
    			}
    		}
    
    		var ele = opt.target;
    		if(typeof opt.target == 'string') ele = document.getElementById(opt.target);
    		var ths = this;
    		shortcut_combination = shortcut_combination.toLowerCase();
    
    		//The function to be called at keypress
    		var func = function(e) {
    			e = e || window.event;
    			
    			if(opt['disable_in_input']) { //Don't enable shortcut keys in Input, Textarea fields
    				var element;
    				if(e.target) element=e.target;
    				else if(e.srcElement) element=e.srcElement;
    				if(element.nodeType==3) element=element.parentNode;
    
    				if(element.tagName == 'INPUT' || element.tagName == 'TEXTAREA') return;
    			}
    	
    			//Find Which key is pressed
    			if (e.keyCode) code = e.keyCode;
    			else if (e.which) code = e.which;
    			var character = String.fromCharCode(code).toLowerCase();
    			
    			if(code == 188) character=","; //If the user presses , when the type is onkeydown
    			if(code == 190) character="."; //If the user presses , when the type is onkeydown
    
    			var keys = shortcut_combination.split("+");
    			//Key Pressed - counts the number of valid keypresses - if it is same as the number of keys, the shortcut function is invoked
    			var kp = 0;
    			
    			//Work around for stupid Shift key bug created by using lowercase - as a result the shift+num combination was broken
    			var shift_nums = {
    				"`":"~",
    				"1":"!",
    				"2":"@",
    				"3":"#",
    				"4":"$",
    				"5":"%",
    				"6":"^",
    				"7":"&",
    				"8":"*",
    				"9":"(",
    				"0":")",
    				"-":"_",
    				"=":"+",
    				";":":",
    				"'":"\"",
    				",":"<",
    				".":">",
    				"/":"?",
    				"\\":"|"
    			}
    			//Special Keys - and their codes
    			var special_keys = {
    				'esc':27,
    				'escape':27,
    				'tab':9,
    				'space':32,
    				'return':13,
    				'enter':13,
    				'backspace':8,
    	
    				'scrolllock':145,
    				'scroll_lock':145,
    				'scroll':145,
    				'capslock':20,
    				'caps_lock':20,
    				'caps':20,
    				'numlock':144,
    				'num_lock':144,
    				'num':144,
    				
    				'pause':19,
    				'break':19,
    				
    				'insert':45,
    				'home':36,
    				'delete':46,
    				'end':35,
    				
    				'pageup':33,
    				'page_up':33,
    				'pu':33,
    	
    				'pagedown':34,
    				'page_down':34,
    				'pd':34,
    	
    				'left':37,
    				'up':38,
    				'right':39,
    				'down':40,
    	
    				'f1':112,
    				'f2':113,
    				'f3':114,
    				'f4':115,
    				'f5':116,
    				'f6':117,
    				'f7':118,
    				'f8':119,
    				'f9':120,
    				'f10':121,
    				'f11':122,
    				'f12':123
    			}
    	
    			var modifiers = { 
    				shift: { wanted:false, pressed:false},
    				ctrl : { wanted:false, pressed:false},
    				alt  : { wanted:false, pressed:false},
    				meta : { wanted:false, pressed:false}	//Meta is Mac specific
    			};
                            
    			if(e.ctrlKey)	modifiers.ctrl.pressed = true;
    			if(e.shiftKey)	modifiers.shift.pressed = true;
    			if(e.altKey)	modifiers.alt.pressed = true;
    			if(e.metaKey)   modifiers.meta.pressed = true;
                            
    			for(var i=0; k=keys[i],i<keys.length; i++) {
    				//Modifiers
    				if(k == 'ctrl' || k == 'control') {
    					kp++;
    					modifiers.ctrl.wanted = true;
    
    				} else if(k == 'shift') {
    					kp++;
    					modifiers.shift.wanted = true;
    
    				} else if(k == 'alt') {
    					kp++;
    					modifiers.alt.wanted = true;
    				} else if(k == 'meta') {
    					kp++;
    					modifiers.meta.wanted = true;
    				} else if(k.length > 1) { //If it is a special key
    					if(special_keys[k] == code) kp++;
    					
    				} else if(opt['keycode']) {
    					if(opt['keycode'] == code) kp++;
    
    				} else { //The special keys did not match
    					if(character == k) kp++;
    					else {
    						if(shift_nums[character] && e.shiftKey) { //Stupid Shift key bug created by using lowercase
    							character = shift_nums[character]; 
    							if(character == k) kp++;
    						}
    					}
    				}
    			}
    			
    			if(kp == keys.length && 
    						modifiers.ctrl.pressed == modifiers.ctrl.wanted &&
    						modifiers.shift.pressed == modifiers.shift.wanted &&
    						modifiers.alt.pressed == modifiers.alt.wanted &&
    						modifiers.meta.pressed == modifiers.meta.wanted) {
    				callback(e);
    	
    				if(!opt['propagate']) { //Stop the event
    					//e.cancelBubble is supported by IE - this will kill the bubbling process.
    					e.cancelBubble = true;
    					e.returnValue = false;
    	
    					//e.stopPropagation works in Firefox.
    					if (e.stopPropagation) {
    						e.stopPropagation();
    						e.preventDefault();
    					}
    					return false;
    				}
    			}
    		}
    		this.all_shortcuts[shortcut_combination] = {
    			'callback':func, 
    			'target':ele, 
    			'event': opt['type']
    		};
    		//Attach the function with the event
    		if(ele.addEventListener) ele.addEventListener(opt['type'], func, false);
    		else if(ele.attachEvent) ele.attachEvent('on'+opt['type'], func);
    		else ele['on'+opt['type']] = func;
    	},
    
    	//Remove the shortcut - just specify the shortcut and I will remove the binding
    	'remove':function(shortcut_combination) {
    		shortcut_combination = shortcut_combination.toLowerCase();
    		var binding = this.all_shortcuts[shortcut_combination];
    		delete(this.all_shortcuts[shortcut_combination])
    		if(!binding) return;
    		var type = binding['event'];
    		var ele = binding['target'];
    		var callback = binding['callback'];
    
    		if(ele.detachEvent) ele.detachEvent('on'+type, callback);
    		else if(ele.removeEventListener) ele.removeEventListener(type, callback, false);
    		else ele['on'+type] = false;
    	}
    }
    
    shortcut.add("Alt+Left",function() {
    	var event = new Event('click');
    	if (settingsShown){
    		document.getElementById('hideSettingsButton_7').dispatchEvent(event);
    	} else{
    		document.getElementById('showSettingsButton_7').dispatchEvent(event);
    	}
    });
    
    if (foldable) moveToBaryCenter();
    
    
    });});
    // COMMON_CODE_BLOCK_END
    </script>
    
    </body>
    </html>



.. figure:: attachment:c3t_maxface.gif
   :alt: {{ :tutorial:c3t_maxface.gif?300 \|}}

   {{ :tutorial:c3t_maxface.gif?300 \|}}

Computing Facets
----------------

An important action that is often needed to come up with new facet
describing inequalities for combinatorial optimization problems is the
computation of convex hulls for small examples.

Pure Integer Case
~~~~~~~~~~~~~~~~~

We begin with the case in which all variables are required to be
integral, i.e., the *pure integer case*. Moreover, the approach depends
on whether the polyhedron is bounded or not.

Bounded Polyhedra
^^^^^^^^^^^^^^^^^

Let us illustrate the approach via the example of the *stable set
problem*: Here one is given an (undirected) Graph G = (V,E) with node
set V and edges E. The goal is to find a largest subset of nodes V’ such
that any two nodes in V’ are not connected by an edge.

For our example consider the 5-cycle, i.e., the graph C5 with five nodes
{1, 2, 3, 4, 5} and edges {1,2}, {2,3}, {3,4}, {4,5}, {5,1}. A
formulation of the stable set problem for this graph looks as follows:

::

   Maximize
    obj: x#1 + x#2 + x#3 + x#4 + x#5
   Subject to
    edge_1: x#2 + x#1 <= 1
    edge_2: x#3 + x#2 <= 1
    edge_3: x#4 + x#3 <= 1
    edge_4: x#5 + x#4 <= 1
    edge_5: x#1 + x#5 <= 1
   Bounds
    0 <= x#1 <= 1
    0 <= x#2 <= 1
    0 <= x#3 <= 1
    0 <= x#4 <= 1
    0 <= x#5 <= 1
   General
    x#1 x#2 x#3 x#4 x#5
   End

Here, ``General`` encodes that the following variables should be
restricted to obtain integer values. Thus, all variables
``x#1, x#2, x#3, x#4, x#5`` are restricted to binary values (0 or 1).
The value 1 encodes that the corresponding node will be in an optimal
stable set (and 0 otherwise). The constraints ``edge_*`` encode that a
most one of the two nodes covered by an edge is selected. Clearly, the
corresponding polyhedron is bounded.

We assume that the above information is contained in the file
``stab.lp``. We now read it into polymake and convert it to rational
form, as explained above:


::

    polymake> $f=lp2poly('files/optimization/stab.lp');
    ........> $p = new Polytope<Rational>($f);

We are now interested in all feasible solutions to the above problem,
i.e., all assignments of 0 or 1 to the variables such that the above
inequalities are satisfied. These feasible points can be computed via:


::

    polymake> $p->LATTICE_POINTS;

To understand these points and make computational use of this
information, we are interested in the convex hull of all feasible
solutions. To this end, we construct a new polytope which is specified
via points for which the convex hull is taken. This can be done as
follows:


::

    polymake> $s=new Polytope(POINTS=>$p->LATTICE_POINTS, COORDINATE_LABELS=>$p->COORDINATE_LABELS);

Here, the coordinate labels, i.e., the variable names, are copied to the
new polytope.

Now, the facets of the new polytope can be computed and listed via:


::

    polymake> print_constraints($s);
    Facets:
    0: x#3 >= 0
    1: -x#4 - x#5 >= -1
    2: x#1 >= 0
    3: -x#2 - x#3 >= -1
    4: -x#1 - x#2 - x#3 - x#4 - x#5 >= -2
    5: -x#1 - x#2 >= -1
    6: x#4 >= 0
    7: -x#1 - x#5 >= -1
    8: -x#3 - x#4 >= -1
    9: x#2 >= 0
    10: x#5 >= 0
    





The facet defining inequalities can be interpreted as follows:

-  There are five trivial inequalities ``x#? >= 0``.

-  The five original ‘edge’ inequalities ``x#i + x#j <= 1`` define
   facets.

-  We have the so-called *odd-cycle inequality*
   ``x#1 +x#2 + x#3 + x#4 + x#5 <= 2``, stating that at most two nodes
   in an (odd) cycle of length 5 can be selected. This inequality can be
   generalized by taking the sum of all variables in an odd cycle and
   restricting the sum to be less or equal to the size of the cycle
   minus 1 divided by 2.

Of course, one can also use the usual polymake output, e.g.,
``print $s->FACETS``.

This example showed one of the routine actions often performed by
discrete optimizers. Of course, this action can also be performed by a
script, which makes the computation a one-line command.

Note that the size of instances that can be handled will probably be
small. Usually, things become difficult from dimension 15 on, but it
depends on the particular structure of your instances, i.e., on the
number of facets and lattice points.

Unbounded Polyhedra
^^^^^^^^^^^^^^^^^^^

If the underlying polyhedron is unbounded, the approach above does not
work anymore, since there are infinitely many lattice points. Arguably,
this case occurs less often than the bounded case, but it is a excellent
show-case for polymake’s potential.

The following mathematical insights are important to treat the unbounded
case. First, we have to assume that the data, i.e., the inequality
description of the polyhedron ``P``, is rational; otherwise, we cannot
expect a finite description of the convex hull. Second, we write ``P``
as the sum of a finite part ``Q`` and the recession cone ``C``. If the
data is rational, the recession cone of ``P`` and of the integer hull
coincide. Third, it suffices to generate the integer points in
``Q + R``, where ``R`` is the parallelotope generated by the rays of
``C``. Thus, ``R`` is generated by the Minkowski sum of the interval
``[0,1]`` and the generating rays.

To illustrate the construction, consider the following example:

::

   Minimize
    obj:  x1 + x2
   Subject to
   C1: x1 + x2 >= 0.5
   C2: x1 - 2 x2 <= 1.5
   C3: x2 - 2 x1 <= 1.5
   General
    x1 x2
   End

.. figure:: attachment:ip-unb.gif
   :alt: {{ :tutorial:ip-unb.gif?300 \| Picture of unbounded polyhedron
   (truncated at upper right)}}

   {{ :tutorial:ip-unb.gif?300 \| Picture of unbounded polyhedron
   (truncated at upper right)}}

We now assume that the example is contained in the file ``unbounded.lp``
and proceed as above


::

    polymake> $f = lp2poly('files/optimization/unbounded.lp');
    ........> $pin = new Polytope<Rational>($f);

The visualization in the picture can be generated with ``$pin->VISUAL``.
The lattice points can be shown with ``$pin->VISUAL->LATTICE_COLORED``.

We now extract the rays of the recession cone


::

    polymake> $rays = $pin->VERTICES->minor($pin->FAR_FACE, All);

This command first computes all vertices of the polyhedron (this
includes unbounded vertices); note that is involves a convex hull
computation. The set ``FAR_FACE`` contains the indices of all vertices
that are unbounded. The result is:


::

    polymake> print $rays;
    0 1 1/2
    0 1 2





Thus, there are two rays that are generators of the recession cone.

We now have to construct the Minkowski hull of all intervals ``[0,r]``
for each ray ``r`` (scaled to be integral). This can be done with the
following code (possibly easier):


::

    polymake> $zero = unit_vector<Rational>($pin->DIM + 1, 0);
    ........> $B = new Polytope<Rational>(POINTS=>$zero);
    ........> foreach my $r (@$rays)
    ........> {
    ........>     $M = new Matrix<Rational>(primitive($r));
    ........>     $M->[0]->[0] = 1;
    ........>     $M = $M / $zero;
    ........>     $ptemp = new Polytope<Rational>(POINTS=>$M);
    ........>     $B = minkowski_sum($B, $ptemp);
    ........> }

The code first generates a polytope ``B`` consisting of 0 only. It then
takes each ray ``r`` in turn and creates ``[0,r]``. It then takes the
Minkowski sum of this new polytope with ``B`` and stores the result in
``B``.

The next step is to obtain the bounded part ``Q`` of ``P``, by first
extracting the bounded vertices and the creating a new polytope:


::

    polymake> $Qpoints = $pin->VERTICES->minor($pin->BOUNDED_VERTICES, All);
    ........> $Q = new Polytope<Rational>(POINTS=>$Qpoints);

The two polytopes are now combined:


::

    polymake> $p = minkowski_sum($Q, $B);

We now generate the lattice points (as in the bounded part) and add the
rays from above:


::

    polymake> $latticemat = new Matrix<Rational>($p->LATTICE_POINTS);
    ........> $newpoints = new Matrix<Rational>($latticemat / $rays);

Here, ``newpoints`` is a matrix that contains all lattice points in
``Q`` and the rays from above.

Finally, the polytope we are interested in is:


::

    polymake> $q = new Polytope(POINTS=>$newpoints, COORDINATE_LABELS=>$pin->COORDINATE_LABELS);

The facets can be viewed as usual:


::

    polymake> print_constraints($q);
    Facets:
    0: 2 x1 - x2 >= -1
    1: 0 >= -1
    2: -x1 + 2 x2 >= -1
    3: x1 + x2 >= 1
    





.. figure:: attachment:ip-unb-integerhull.gif
   :alt: {{ :tutorial:ip-unb-integerhull.gif?300 \|}}

   {{ :tutorial:ip-unb-integerhull.gif?300 \|}}

Note that the upper right part (including the red vertices) arises from
truncation of the polyhedron for visualization.

Mixed-Integer Case
~~~~~~~~~~~~~~~~~~

Let us now briefly discuss how to proceed if there are variables that
are allowed to be integral. In this case there are several different
types of information that one might be interested in. Let us first
consider the question of how to compute the convex hull of all feasible
integral variables, i.e., we consider the projection to the integral
variables and then consider the convex hull of all feasible solutions.
We only consider the bounded case, i.e., the original polyhedron is
bounded.

Consider the following example:

::

   Minimize
    obj:  x1 + x2
   Subject to
   C1: s1 - 10 x1 <= 0
   C2: s2 - 10 x2 <= 0
   C3: s1 + s2 <= 1.5
   C4: s1 + s2 >= 0.5
   Bounds
    0 <= s1
    0 <= s2
    0 <= x1 <= 1
    0 <= x2 <= 1
   General
    x1 x2
   End

In this example there are two integral variables ``x1`` and ``x2``,
while ``s1`` and ``s2`` are continuous variables. Assuming the data is
contained in the file ``mip.lp``, we proceed as follows:


::

    polymake> $m=lp2poly('files/optimization/mip.lp');
    ........> $p = new Polytope<Rational>($m);

We project the polyhedron in ``$p`` to the third and fourth variables as
follows:


::

    polymake> $q=projection($p, [3,4]);

We now construct the convex hull of all feasible points as above:


::

    polymake> $s=new Polytope(POINTS=>$q->LATTICE_POINTS);
    ........> print_constraints($s);
    Facets:
    0: -x1 >= -1
    1: -x2 >= -1
    2: x1 + x2 >= 1
    





Thus, as expected, the convex hull equals the triangle with vertices
``{(0,1),(1,0),(1,1)}``.

Integral Polytopes and Total Unimodularity
------------------------------------------

As explained in the previous example, the integral points in a polytope
are of particular interest in discrete optimization. These points are
called *lattice points* in polymake and the corresponding convex hull
*lattice polytope*. The handling of such polytopes is explained in more
detail in the `Tutorial for Lattice
Polytopes <lattice_polytopes_tutorial>`__.

Of particular interest for discrete optimization are properties of the
original inequality system to define a lattice polytope, i.e., a
polytope such that all of its vertices are integral (this can be tested
by checking the property ``LATTICE``). One particularly interesting case
occurs if the matrix defining the polytope is *totally unimodular* and
the right hand side is integral.

Using the polymake extension
`Unimodularity <https///github.com/xammy/unimodularity-test/wiki/Polymake-Extension>`__
by Matthias Walter, this can be checked as illustrated in the following
examples.

Example: Explicit Matrix
~~~~~~~~~~~~~~~~~~~~~~~~

In a first example, we directly create an integral matrix


::

    polymake> $M=new Matrix<Integer>([[1,1,0,0],[1,0,1,0],[1,0,0,1]]);

The total unimodularity of this matrix can be checked as follows:


::

    polymake> print totally_unimodular($M);
    true




Thus, the given matrix is totally unimodular.

Example: Matrix from Input
~~~~~~~~~~~~~~~~~~~~~~~~~~

In the second example, we reuse the file ``c3t`` from the example above.
We read it into polymake:


::

    polymake> $f=lp2poly('files/optimization/c3t.lp');
    ........> $p = new Polytope<Rational>($f);

We now want to check whether the constraint matrix defined by the
inequalities is totally unimodular (note that there are no equations in
this example). Thus we first extract the inequality matrix without the
first column (as an integer matrix) and then perform the test:


::

    polymake> $A = new Matrix<Integer>($p->INEQUALITIES->minor(All, ~[0]));
    ........> print totally_unimodular($A);
    true




Thus, this matrix is totally unimodular as well.

Total dual integrality
----------------------

Computations with respect to total dual integrality (TDI) can also be
performed in polymake. Currently (August 2013), you need the perpetual
beta version of polymake to access this functionality.

The main functions are:

-  The function ``totally_dual_integral`` takes an inequality system (as
   a matrix) and checks whether it is totally dual integral.

-  The function ``make_totally_dual_integral`` takes a polytope and
   returns a new polytope with inequalities that are TDI.

Note that the input has to be full-dimensional in order to use these
functions.

To demonstrate the behavior of these functions, consider the 5-cycle
example from above again:

::

   Maximize
    obj: x#1 + x#2 + x#3 + x#4 + x#5
   Subject to
    edge_1: x#2 + x#1 <= 1
    edge_2: x#3 + x#2 <= 1
    edge_3: x#4 + x#3 <= 1
    edge_4: x#5 + x#4 <= 1
    edge_5: x#1 + x#5 <= 1
   Bounds
    0 <= x#1 <= 1
    0 <= x#2 <= 1
    0 <= x#3 <= 1
    0 <= x#4 <= 1
    0 <= x#5 <= 1
   General
    x#1 x#2 x#3 x#4 x#5
   End

Let us test whether the inequality system of this example is TDI. Thus,
we first load the data as usual:


::

    polymake> $f = lp2poly('files/optimization/stab.lp');
    ........> $p = new Polytope<Rational>($f);

We now extract the corresponding inequality system and check it for
TDIness:


::

    polymake> $M = new Matrix<Rational>($p->INEQUALITIES);
    ........> print totally_dual_integral($M);
    false




.. raw:: html

    <details><summary><pre style="display:inline"><small>Click here for additional output</small></pre></summary>
    <pre>
    polymake: used package libnormaliz
      Normaliz is a tool for computations in affine monoids, vector configurations, lattice polytopes, and rational cones.
      Copyright by Winfried Bruns, Bogdan Ichim, Christof Soeger.
      http://www.math.uos.de/normaliz/
    
    </pre>
    </details>




The system is not TDI, which we expected from general theory, since we
know that the polytope is not integral, but the system has integral
coefficients. Consequently, let us construct a TDI-system for this
polytope:


::

    polymake> $q = make_totally_dual_integral($p);
    ........> print_constraints($q);
    Inequalities:
    0: x5 >= 0
    1: x4 >= 0
    2: x3 >= 0
    3: x2 >= 0
    4: x1 >= 0
    5: -x1 - x2 >= -1
    6: -x1 - x5 >= -1
    7: -x2 - x3 >= -1
    8: -x3 - x4 >= -1
    9: -x4 - x5 >= -1
    10: -x1 - x2 - x3 - x4 - x5 >= -5/2
    11: 0 >= -1
    





As expected, the right hand side is non integral (otherwise, we know
from general theory that the polytope would be integral as well). The
result is now TDI:


::

    polymake> $N = new Matrix<Rational>($q->INEQUALITIES);
    ........> print totally_dual_integral($N);
    true




*Note* that we need to take the inequalities instead of facets here,
since facets are irredundant and thus might not be TDI, although the
complete set of inequalities is TDI.

Chvátal-Gomory Closure
----------------------

In the following we want to briefly show how closures of polytopes with
respect to certain cutting-plane operators can be computed. We consider
the two well-known cutting-plane operators here. The first one is the
Chvátal-Gomory generator and the second one is the Lift-and-project
operator as defined by Balas. For simplicity we will assume that the
considered polytope is full-dimensional.

Chvátal-Gomory Closure - Example 1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We first consider the polytope from the stable set problem from above:

::

   Maximize
    obj: x#1 + x#2 + x#3 + x#4 + x#5
   Subject to
    edge_1: x#2 + x#1 <= 1
    edge_2: x#3 + x#2 <= 1
    edge_3: x#4 + x#3 <= 1
    edge_4: x#5 + x#4 <= 1
    edge_5: x#1 + x#5 <= 1
   Bounds
    0 <= x#1 <= 1
    0 <= x#2 <= 1
    0 <= x#3 <= 1
    0 <= x#4 <= 1
    0 <= x#5 <= 1
   General
    x#1 x#2 x#3 x#4 x#5
   End

As before we read in the file using ``lp2poly``:


::

    polymake> $f = lp2poly('files/optimization/stab.lp');
    ........> $p = new Polytope<Rational>($f);

The Chvátal-Gomory closure of a polytope can be computed with the
function ``gc_closure``. The function takes a full-dimensional polytope
and returns a new polytope. This contains the system of inequalities
defining the closure in the property ``INEQUALITIES``. For our example,
we obtain:


::

    polymake> $g = gc_closure($p);
    ........> print print_constraints($g);
    Inequalities:
    0: x5 >= 0
    1: x4 >= 0
    2: x3 >= 0
    3: x2 >= 0
    4: x1 >= 0
    5: -x1 - x2 >= -1
    6: -x1 - x5 >= -1
    7: -x2 - x3 >= -1
    8: -x3 - x4 >= -1
    9: -x4 - x5 >= -1
    10: -x1 - x2 - x3 - x4 - x5 >= -2
    11: 0 >= -1
    





Let us check whether the resulting polytope is integral:


::

    polymake> print $g->LATTICE;
    true




Thus, in this case, we have obtained the integer hull by one step of the
Chvatal-Gomory-closure.

Chvátal-Gomory Closure - Example 2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us now consider the classical example of a polytope with the
vertices of simplex in d dimensions and the point 1/2 times (1, …, 1).
It can be shown that such a polytope has rank at least log(d) - 1, see
`Pokutta, 2011 <http://www.box.net/shared/at1y8i3pq434bxt6m9xm>`__]. In
our example, we use d = 4:


::

    polymake> $M = new Matrix<Rational>([[1,0,0,0,0],[1,1,0,0,0],[1,0,1,0,0],[1,0,0,1,0],[1,0,0,0,1],[1,1/2,1/2,1/2,1/2]]);
    ........> $t = new Polytope<Rational>(POINTS => $M);
    ........> $t1 = gc_closure($t);
    ........> $t1->FACETS;
    ........> print_constraints($t1);
    Facets:
    0: x4 >= 0
    1: x3 >= 0
    2: x2 >= 0
    3: x1 >= 0
    4: -x1 - x2 - x3 >= -1
    5: -x1 - x2 - x4 >= -1
    6: -x1 - x3 - x4 >= -1
    7: -x2 - x3 - x4 >= -1
    





::

    polymake> print $t1->LATTICE;
    false




Thus, one round was not enough to produce an integral polytope. Indeed,
the vertices are


::

    polymake> $t1->VERTICES;
    ........> print $t1->VERTICES;
    1 1 0 0 0
    1 0 0 0 0
    1 1/3 1/3 1/3 1/3
    1 0 1 0 0
    1 0 0 1 0
    1 0 0 0 1





However, one more round is enough:


::

    polymake> $t2 = gc_closure($t1);
    ........> $t2->FACETS;
    ........> print_constraints($t2);
    Facets:
    0: x4 >= 0
    1: x3 >= 0
    2: x2 >= 0
    3: x1 >= 0
    4: -x1 - x2 - x3 - x4 >= -1
    





::

    polymake> print $t2->LATTICE;
    true




Lift-and-project closure
------------------------

The lift-and-project closure of a 0/1-polytope P can be generated as
follows: for each coordinate compute the intersection of P with the pair
of opposite cube faces and compute the convex hull. Then intersect the
result with P. The following script defines a subroutine performing this
operation - the code is somewhat complicated throught the fact that we
need to check whether parts are empty.


::

   use application "polytope";

   sub lpclosure
   {
       my $p = shift;
       my $d = $p->AMBIENT_DIM;
       my $q = new Polytope<Rational>($p);
       for (my $k = 0; $k < $d; $k = $k+1)
       {
           if ( $q->DIM == -1 )         # can stop as soon as $q is empty
           {
                return $q;
           }
       
           # create reversed opposite inequalities of 0/1-cube and corresponding polyhedra
           my $v1 = new Vector<Rational>(0 | -unit_vector($d, $k));
           my $v2 = new Vector<Rational>(-1 | unit_vector($d, $k));
       
           # create intersection of corresponding polyhedra with iterated polyhedron $q
           my $b1 = new Polytope<Rational>(INEQUALITIES => $v1 / $q->FACETS);
           my $b2 = new Polytope<Rational>(INEQUALITIES => $v2 / $q->FACETS);
       
           if ( ($b1->DIM > -1) && ($b2->DIM > -1) )
           {
               my $c = conv($b1, $b2);
               $q = intersection($q, $c);
           }
           elsif ( ($b1->DIM > -1) && ($b2->DIM == -1) )
           {
               $q = intersection($q, $b1);
           }
           elsif ( ($b1->DIM == -1) && ($b2->DIM > -1) )
           {
               $q = intersection($q, $b2);
           }
       }
       return $q;
   }


::

    polymake> script("files/optimization/lpclosure.pl");

Lift-and-Project Closure - Example 1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For our well known stable set example, we get the following:


::

    polymake> $q = lpclosure($p);
    ........> $q->FACETS;
    ........> print_constraints($q);
    Facets:
    0: -x1 - x2 >= -1
    1: x3 >= 0
    2: x2 >= 0
    3: x4 >= 0
    4: -x3 - x4 >= -1
    5: -x4 - x5 >= -1
    6: x1 >= 0
    7: -x2 - x3 >= -1
    8: -x1 - x2 - x3 - x4 - x5 >= -2
    9: -x1 - x5 >= -1
    10: x5 >= 0
    





Thus, the lift-and-project closure in this case gives the integral hull
(as we have seen above).

Lift-and-Project Closure - Example 2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us now consider the same example as for CG-closures:


::

    polymake> $M = new Matrix<Rational>([[1,0,0,0,0],[1,1,0,0,0],[1,0,1,0,0],[1,0,0,1,0],[1,0,0,0,1],[1,1/2,1/2,1/2,1/2]]);
    ........> $p = new Polytope<Rational>(POINTS => $M);
    ........> $q = lpclosure($p);
    ........> $q->FACETS;
    ........> print_constraints($q);
    Facets:
    0: x1 >= 0
    1: x4 >= 0
    2: -x1 - x2 - x3 - x4 >= -1
    3: x2 >= 0
    4: x3 >= 0
    





Thus, we have obtained the integral hull in a single step of the
lift-and-project closure as opposed to two steps in the CG-closure.
