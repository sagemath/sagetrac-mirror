.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Regular subdivisons
-------------------

*Regular subdivsions* of point sets appear in several different
applications. ``polymake`` allows to define regular subdivisions of
point configurations (e.g. the lattice points of a lattice polytope) via
weights on the points. The weights define a *height function* on the
points, and the subdivision is described by the lower hull of the
polytope defined by the lifted points (see
`here <http://link.springer.com/book/10.1007%2F978-3-642-12971-1>`__ for
example to get an idea of the mathematical backround). Let us look at an
example on how to create a regular subdivision.


.. link

.. CODE-BLOCK:: perl

    polymake> $M = new Matrix<Rational>([[1,0,0],[1,2,0],[1,1,1],[1,0,2],[1,2,2],[1,1,3]]);
    polymake> $w = new Vector<Rational>([0,0,1,0,0,3]);
    polymake> $S = new fan::SubdivisionOfPoints(POINTS=>$M,WEIGHTS=>$w);
    polymake> print $S->MAXIMAL_CELLS;
    {3 4 5}
    {0 1 3 4}





.. raw:: html

    <details><summary><pre style="display:inline"><small>Click here for additional output</small></pre></summary>
    <pre>
    polymake: used package cdd
      cddlib
      Implementation of the double description method of Motzkin et al.
      Copyright by Komei Fukuda.
      http://www-oldurls.inf.ethz.ch/personal/fukudak/cdd_home/
    
    </pre>
    </details>




::

   {0 1 3 4}
   {3 4 5}

The six points (given by the six rows of the matrix ``M``) define a
pentagon with interior point ``(1,1)``. The weights ``w`` must be given
in the same order as the vertices. So in our case the interior point has
weight ``1``. Switching to the ``fan`` application we may create the
regular subdivision via the polymake object ``SubdivisionOfPoints``,
whose input parameters are the points ``M`` and weights ``w``. A maximal
cell of the subdivision is given by a set of indices, representing the
points that cell contains. In our example we got two maximal cells and
the first one is a quadrilateral with vertices 0, 1, 3 and 4.

We may visualize the regular subdivision if it is at most
three-dimensional.


.. link

.. CODE-BLOCK:: perl

    polymake> $S->VISUAL;


.. raw:: html

    <!--
    polymake for andrew
    Thu Mar 28 11:50:05 2019
    pcom:S
    -->
    
    
    <html>
       <head>
          <title>pcom:S</title>
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
    
    		<div id='settings_0' class='settings'>
    			<div class=group id='explode_0'>
    				<strong>Explode</strong>
    				<input id='explodeRange_0' type='range' min=0 max=6 step=0.01 value=0>
    				<div class=indented><input id='explodeCheckbox_0' type='checkbox'>Automatic explosion</div>
    				<div class=suboption>Exploding speed</div>
    				<input id='explodingSpeedRange_0' type='range' min=0 max=0.5 step=0.001 value=0.05>
    			</div>
    
    			
    			<div class=group id='transparency_0' class='transparency'>
    				<strong>Transparency</strong>
    				<input id='transparencyRange_0' type='range' min=0 max=1 step=0.01 value=0>
    			</div>
    			
    			<div class=group id='rotation_0'>
    				<strong>Rotation</strong>
    				<div class=indented>
    					<div><input type='checkbox' id='changeRotationX_0'> x-axis</div>
    					<div><input type='checkbox' id='changeRotationY_0'> y-axis</div>
    					<div><input type='checkbox' id='changeRotationZ_0'> z-axis</div>
    					<button id='resetButton_0' class='resetButton' >Reset</button>
    				</div>
    
    				<div class=suboption>Rotation speed</div>
    				<input id='rotationSpeedRange_0' type='range' min=0 max=5 step=0.01 value=2>
    
    			</div>
    
    
    			<div class=group id='display_0'>
    				<strong>Display</strong>
    				<div class=indented>
    					<div id='shownObjectsList_0' class='shownObjectsList'></div>
    					<div class='labelsCheckbox'><input type='checkbox' id='labelsCheckboxInput_0' checked>Labels</div>
    				</div>
    			</div>
    
    
    			<div class=group id='svg_0'>
    				<strong>SVG</strong>
    				<div class=indented>
    					<form>
    						<input type="radio" name='screenshotMode' value='download' id='download_0' checked> Download<br>
    						<input type="radio" name='screenshotMode' value='tab' id='tab_0' > New tab<br>
    					</form>
    					<button id='takeScreenshot_0'>Screenshot</button>
    				</div>
    			</div>
    
    		</div>	<!-- end of settings -->
    		<img id='hideSettingsButton_0' style="display: none" class='hideSettingsButton' src='/kernelspecs/polymake/close.svg' width=20px">
    		<img id='showSettingsButton_0' class='showSettingsButton' src='/kernelspecs/polymake/menu.svg' width=20px">
    <div id="model85250885745"></div>
    
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
       var container = document.getElementById( 'model85250885745' );
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
       var target = document.querySelector('#model85250885745');
    
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
    
       var objectnames = ["points","Vertices","unnamed__1","unnamed__2"];
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 0, 0));
       allpoints.push(new PMPoint(2, 0, 0));
       allpoints.push(new PMPoint(1, 1, 0));
       allpoints.push(new PMPoint(0, 2, 0));
       allpoints.push(new PMPoint(2, 2, 0));
       allpoints.push(new PMPoint(1, 3, 0));
    
       <!-- Vertex style -->
       var points_material = new THREE.MeshBasicMaterial ( {color: 0x000000, } );
    
       points_material.side = THREE.DoubleSide;
       points_material.transparent = true;
    
       <!-- POINTS -->
       allpoints[0].makesphere(0.02,points_material);
       allpoints[1].makesphere(0.02,points_material);
       allpoints[2].makesphere(0.02,points_material);
       allpoints[3].makesphere(0.02,points_material);
       allpoints[4].makesphere(0.02,points_material);
       allpoints[5].makesphere(0.02,points_material);
       allpoints[0].makelabel("0");
       allpoints[1].makelabel("1");
       allpoints[2].makelabel("2");
       allpoints[3].makelabel("3");
       allpoints[4].makelabel("4");
       allpoints[5].makelabel("5");
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 0, 0));
       allpoints.push(new PMPoint(2, 0, 0));
       allpoints.push(new PMPoint(0, 2, 0));
       allpoints.push(new PMPoint(2, 2, 0));
       allpoints.push(new PMPoint(1, 3, 0));
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       <!-- Edge style -->
       var line_material = new THREE.LineBasicMaterial ( {color: 0x000000, linewidth: 1.5, } );
    
       line_material.side = THREE.DoubleSide;
       line_material.transparent = true;
    
       <!-- EDGES -->
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[3].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 2, 0));
       allpoints.push(new PMPoint(2, 2, 0));
       allpoints.push(new PMPoint(1, 3, 0));
    
       <!-- Vertex style -->
       var materials = [
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
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 0.5, side: THREE.DoubleSide , depthWrite: false, depthTest: false, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(2, 0, 1, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 0, 0));
       allpoints.push(new PMPoint(2, 0, 0));
       allpoints.push(new PMPoint(0, 2, 0));
       allpoints.push(new PMPoint(2, 2, 0));
    
       <!-- Vertex style -->
       var materials = [
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
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
       faces.vertices.push(allpoints[3].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 0.5, side: THREE.DoubleSide , depthWrite: false, depthTest: false, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(1, 3, 2, undefined, undefined, 0));
       faces.faces.push(new THREE.Face3(1, 2, 0, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
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
    		document.getElementById('settings_0').style.position = 'absolute';
    		document.getElementById('settings_0').style.display = 'block';
    		document.getElementById('showSettingsButton_0').style.display = 'none';
    		document.getElementById('hideSettingsButton_0').style.display = 'block';
    		settingsShown = true;
    	}
    
    	function hideSettings(event){
    		event.currentTarget.style.display = 'none';
    		document.getElementById('settings_0').style.display = 'none';
    		document.getElementById('hideSettingsButton_0').style.display = 'none';
    		document.getElementById('showSettingsButton_0').style.display = 'block';
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
    
    		if (document.getElementById('tab_0').checked){
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
    							document.getElementById('transparency_0').style.display = 'block';
    							document.getElementById('transparencyRange_0').value = 1 - 
    								child.material.materials[k].opacity;
    							return;
    						}
    					}
    				} else if (	child.material.transparent && 
    								child.material.type == "MeshBasicMaterial" &&
    								child.geometry.type == "Geometry"){
    					document.getElementById('transparency_0').style.display = 'block';
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
    		document.getElementById('explodeRange_0').value = explodeScale;
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
    		document.getElementById('explodeRange_0').value = explodeScale;
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
    	var shownObjectsList = document.getElementById('shownObjectsList_0');
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
    tempobj = document.getElementById('explodeRange_0');
    if (tempobj) {
       tempobj.oninput = triggerExplode;
       document.getElementById('explodeCheckbox_0').onchange = triggerAutomaticExplode;
       document.getElementById('explodingSpeedRange_0').oninput = setExplodingSpeed;
    }
    tempobj = document.getElementById('foldRange_0');
    if (tempobj) {
       tempobj.oninput = fold;
    }
    document.getElementById('transparencyRange_0').oninput = changeTransparency;
    document.getElementById('changeRotationX_0').onchange = changeRotationX;
    document.getElementById('changeRotationY_0').onchange = changeRotationY;
    document.getElementById('changeRotationZ_0').onchange = changeRotationZ;
    document.getElementById('resetButton_0').onclick = resetScene;
    document.getElementById('rotationSpeedRange_0').oninput = changeRotationSpeedFactor;
    document.getElementById('labelsCheckboxInput_0').onchange = displayOrHideLabels;
    document.getElementById('takeScreenshot_0').onclick = takeSvgScreenshot;
    document.getElementById('showSettingsButton_0').onclick = showSettings;
    document.getElementById('hideSettingsButton_0').onclick = hideSettings;
    
    	
    
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
    		document.getElementById('hideSettingsButton_0').dispatchEvent(event);
    	} else{
    		document.getElementById('showSettingsButton_0').dispatchEvent(event);
    	}
    });
    
    if (foldable) moveToBaryCenter();
    
    
    });});
    // COMMON_CODE_BLOCK_END
    </script>
    
    </body>
    </html>



.. raw:: html

    <details><summary><pre style="display:inline"><small>Click here for additional output</small></pre></summary>
    <pre>
    polymake: used package threejs
       Three.js is a lightweight cross-browser JavaScript library/API used to create and display animated 3D computer graphics on a Web browser.
       See http://github.com/mrdoob for the source code.
    
    </pre>
    </details>




If you use javaview for visualization, then this might look similar to
the following:

.. figure:: attachment:subdiv_of_points_1.gif
   :alt: {{:tutorial:subdiv_of_points_1.gif?300|}}

   {{:tutorial:subdiv_of_points_1.gif?300|}}

Note that the quadrilateral contains point 2 in its interior and that
this point is colored black. This corresponds to the fact that the
lifted point 2 lies above the convex hull of the lifted points 0,1,3 and
4. Therefore the maximal cell describing the quadrilateral does not
contain the point 2. We may change the lifting function by giving point
2 height ``0`` as well:


.. link

.. CODE-BLOCK:: perl

    polymake> $w = new Vector<Rational>([0,0,0,0,0,3]);
    polymake> $S = new fan::SubdivisionOfPoints(POINTS=>$M,WEIGHTS=>$w);
    polymake> print $S->MAXIMAL_CELLS;
    {3 4 5}
    {0 1 2 3 4}





::

   {0 1 2 3 4}
   {3 4 5}


.. link

.. CODE-BLOCK:: perl

    polymake> $S->VISUAL;


.. raw:: html

    <!--
    polymake for andrew
    Thu Mar 28 11:52:14 2019
    pcom:S
    -->
    
    
    <html>
       <head>
          <title>pcom:S</title>
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
    
    		<div id='settings_1' class='settings'>
    			<div class=group id='explode_1'>
    				<strong>Explode</strong>
    				<input id='explodeRange_1' type='range' min=0 max=6 step=0.01 value=0>
    				<div class=indented><input id='explodeCheckbox_1' type='checkbox'>Automatic explosion</div>
    				<div class=suboption>Exploding speed</div>
    				<input id='explodingSpeedRange_1' type='range' min=0 max=0.5 step=0.001 value=0.05>
    			</div>
    
    			
    			<div class=group id='transparency_1' class='transparency'>
    				<strong>Transparency</strong>
    				<input id='transparencyRange_1' type='range' min=0 max=1 step=0.01 value=0>
    			</div>
    			
    			<div class=group id='rotation_1'>
    				<strong>Rotation</strong>
    				<div class=indented>
    					<div><input type='checkbox' id='changeRotationX_1'> x-axis</div>
    					<div><input type='checkbox' id='changeRotationY_1'> y-axis</div>
    					<div><input type='checkbox' id='changeRotationZ_1'> z-axis</div>
    					<button id='resetButton_1' class='resetButton' >Reset</button>
    				</div>
    
    				<div class=suboption>Rotation speed</div>
    				<input id='rotationSpeedRange_1' type='range' min=0 max=5 step=0.01 value=2>
    
    			</div>
    
    
    			<div class=group id='display_1'>
    				<strong>Display</strong>
    				<div class=indented>
    					<div id='shownObjectsList_1' class='shownObjectsList'></div>
    					<div class='labelsCheckbox'><input type='checkbox' id='labelsCheckboxInput_1' checked>Labels</div>
    				</div>
    			</div>
    
    
    			<div class=group id='svg_1'>
    				<strong>SVG</strong>
    				<div class=indented>
    					<form>
    						<input type="radio" name='screenshotMode' value='download' id='download_1' checked> Download<br>
    						<input type="radio" name='screenshotMode' value='tab' id='tab_1' > New tab<br>
    					</form>
    					<button id='takeScreenshot_1'>Screenshot</button>
    				</div>
    			</div>
    
    		</div>	<!-- end of settings -->
    		<img id='hideSettingsButton_1' style="display: none" class='hideSettingsButton' src='/kernelspecs/polymake/close.svg' width=20px">
    		<img id='showSettingsButton_1' class='showSettingsButton' src='/kernelspecs/polymake/menu.svg' width=20px">
    <div id="model10090275001"></div>
    
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
       var container = document.getElementById( 'model10090275001' );
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
       var target = document.querySelector('#model10090275001');
    
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
    
       var objectnames = ["points","Vertices","unnamed__1","unnamed__2"];
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 0, 0));
       allpoints.push(new PMPoint(2, 0, 0));
       allpoints.push(new PMPoint(1, 1, 0));
       allpoints.push(new PMPoint(0, 2, 0));
       allpoints.push(new PMPoint(2, 2, 0));
       allpoints.push(new PMPoint(1, 3, 0));
    
       <!-- Vertex style -->
       var points_material = new THREE.MeshBasicMaterial ( {color: 0x000000, } );
    
       points_material.side = THREE.DoubleSide;
       points_material.transparent = true;
    
       <!-- POINTS -->
       allpoints[0].makesphere(0.02,points_material);
       allpoints[1].makesphere(0.02,points_material);
       allpoints[2].makesphere(0.02,points_material);
       allpoints[3].makesphere(0.02,points_material);
       allpoints[4].makesphere(0.02,points_material);
       allpoints[5].makesphere(0.02,points_material);
       allpoints[0].makelabel("0");
       allpoints[1].makelabel("1");
       allpoints[2].makelabel("2");
       allpoints[3].makelabel("3");
       allpoints[4].makelabel("4");
       allpoints[5].makelabel("5");
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 0, 0));
       allpoints.push(new PMPoint(2, 0, 0));
       allpoints.push(new PMPoint(0, 2, 0));
       allpoints.push(new PMPoint(2, 2, 0));
       allpoints.push(new PMPoint(1, 3, 0));
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       <!-- Edge style -->
       var line_material = new THREE.LineBasicMaterial ( {color: 0x000000, linewidth: 1.5, } );
    
       line_material.side = THREE.DoubleSide;
       line_material.transparent = true;
    
       <!-- EDGES -->
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[3].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 2, 0));
       allpoints.push(new PMPoint(2, 2, 0));
       allpoints.push(new PMPoint(1, 3, 0));
    
       <!-- Vertex style -->
       var materials = [
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
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 0.5, side: THREE.DoubleSide , depthWrite: false, depthTest: false, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(2, 0, 1, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 0, 0));
       allpoints.push(new PMPoint(2, 0, 0));
       allpoints.push(new PMPoint(1, 1, 0));
       allpoints.push(new PMPoint(0, 2, 0));
       allpoints.push(new PMPoint(2, 2, 0));
    
       <!-- Vertex style -->
       var materials = [
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFFFF00, }),
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
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 0.5, side: THREE.DoubleSide , depthWrite: false, depthTest: false, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(1, 4, 3, undefined, undefined, 0));
       faces.faces.push(new THREE.Face3(1, 3, 0, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
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
    		document.getElementById('settings_1').style.position = 'absolute';
    		document.getElementById('settings_1').style.display = 'block';
    		document.getElementById('showSettingsButton_1').style.display = 'none';
    		document.getElementById('hideSettingsButton_1').style.display = 'block';
    		settingsShown = true;
    	}
    
    	function hideSettings(event){
    		event.currentTarget.style.display = 'none';
    		document.getElementById('settings_1').style.display = 'none';
    		document.getElementById('hideSettingsButton_1').style.display = 'none';
    		document.getElementById('showSettingsButton_1').style.display = 'block';
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
    
    		if (document.getElementById('tab_1').checked){
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
    							document.getElementById('transparency_1').style.display = 'block';
    							document.getElementById('transparencyRange_1').value = 1 - 
    								child.material.materials[k].opacity;
    							return;
    						}
    					}
    				} else if (	child.material.transparent && 
    								child.material.type == "MeshBasicMaterial" &&
    								child.geometry.type == "Geometry"){
    					document.getElementById('transparency_1').style.display = 'block';
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
    		document.getElementById('explodeRange_1').value = explodeScale;
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
    		document.getElementById('explodeRange_1').value = explodeScale;
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
    	var shownObjectsList = document.getElementById('shownObjectsList_1');
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
    tempobj = document.getElementById('explodeRange_1');
    if (tempobj) {
       tempobj.oninput = triggerExplode;
       document.getElementById('explodeCheckbox_1').onchange = triggerAutomaticExplode;
       document.getElementById('explodingSpeedRange_1').oninput = setExplodingSpeed;
    }
    tempobj = document.getElementById('foldRange_1');
    if (tempobj) {
       tempobj.oninput = fold;
    }
    document.getElementById('transparencyRange_1').oninput = changeTransparency;
    document.getElementById('changeRotationX_1').onchange = changeRotationX;
    document.getElementById('changeRotationY_1').onchange = changeRotationY;
    document.getElementById('changeRotationZ_1').onchange = changeRotationZ;
    document.getElementById('resetButton_1').onclick = resetScene;
    document.getElementById('rotationSpeedRange_1').oninput = changeRotationSpeedFactor;
    document.getElementById('labelsCheckboxInput_1').onchange = displayOrHideLabels;
    document.getElementById('takeScreenshot_1').onclick = takeSvgScreenshot;
    document.getElementById('showSettingsButton_1').onclick = showSettings;
    document.getElementById('hideSettingsButton_1').onclick = hideSettings;
    
    	
    
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
    		document.getElementById('hideSettingsButton_1').dispatchEvent(event);
    	} else{
    		document.getElementById('showSettingsButton_1').dispatchEvent(event);
    	}
    });
    
    if (foldable) moveToBaryCenter();
    
    
    });});
    // COMMON_CODE_BLOCK_END
    </script>
    
    </body>
    </html>



.. figure:: attachment:subdiv_of_points_2.gif
   :alt: {{:tutorial:subdiv_of_points_2.gif?300|}}

   {{:tutorial:subdiv_of_points_2.gif?300|}}

This new height function generates the same polyhedral complex, yet a
different subdivision of the point set. Note that the maximal cell that
describes the quadrilateral now contains the point 2, but it is not a
vertex of that cell. This is the reason for the yellow coloring. This
corresponds to the fact that its lifting lies in the convex hull of the
lifted quadrilateral, whereat it is not a vertex of the lifted polytope.
Changing the height of point 2 again to ``-1`` the subdivision should
become a triangulation of the point set.


.. link

.. CODE-BLOCK:: perl

    polymake> $w = new Vector<Rational>([0,0,-1,0,0,3]);
    polymake> $S = new fan::SubdivisionOfPoints(POINTS=>$M,WEIGHTS=>$w);
    polymake> print $S->MAXIMAL_CELLS;
    {0 1 2}
    {0 2 3}
    {2 3 4}
    {3 4 5}
    {1 2 4}





::

   {0 2 3}
   {0 1 2}
   {1 2 4}
   {2 3 4}
   {3 4 5}


.. link

.. CODE-BLOCK:: perl

    polymake> $S->VISUAL;


.. raw:: html

    <!--
    polymake for andrew
    Thu Apr  4 10:31:25 2019
    pcom:S
    -->
    
    
    <html>
       <head>
          <title>pcom:S</title>
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
    
    		<div id='settings_0' class='settings'>
    			<div class=group id='explode_0'>
    				<strong>Explode</strong>
    				<input id='explodeRange_0' type='range' min=0 max=6 step=0.01 value=0>
    				<div class=indented><input id='explodeCheckbox_0' type='checkbox'>Automatic explosion</div>
    				<div class=suboption>Exploding speed</div>
    				<input id='explodingSpeedRange_0' type='range' min=0 max=0.5 step=0.001 value=0.05>
    			</div>
    
    			
    			<div class=group id='transparency_0' class='transparency'>
    				<strong>Transparency</strong>
    				<input id='transparencyRange_0' type='range' min=0 max=1 step=0.01 value=0>
    			</div>
    			
    			<div class=group id='rotation_0'>
    				<strong>Rotation</strong>
    				<div class=indented>
    					<div><input type='checkbox' id='changeRotationX_0'> x-axis</div>
    					<div><input type='checkbox' id='changeRotationY_0'> y-axis</div>
    					<div><input type='checkbox' id='changeRotationZ_0'> z-axis</div>
    					<button id='resetButton_0' class='resetButton' >Reset</button>
    				</div>
    
    				<div class=suboption>Rotation speed</div>
    				<input id='rotationSpeedRange_0' type='range' min=0 max=5 step=0.01 value=2>
    
    			</div>
    
    
    			<div class=group id='display_0'>
    				<strong>Display</strong>
    				<div class=indented>
    					<div id='shownObjectsList_0' class='shownObjectsList'></div>
    					<div class='labelsCheckbox'><input type='checkbox' id='labelsCheckboxInput_0' checked>Labels</div>
    				</div>
    			</div>
    
    
    			<div class=group id='svg_0'>
    				<strong>SVG</strong>
    				<div class=indented>
    					<form>
    						<input type="radio" name='screenshotMode' value='download' id='download_0' checked> Download<br>
    						<input type="radio" name='screenshotMode' value='tab' id='tab_0' > New tab<br>
    					</form>
    					<button id='takeScreenshot_0'>Screenshot</button>
    				</div>
    			</div>
    
    		</div>	<!-- end of settings -->
    		<img id='hideSettingsButton_0' style="display: none" class='hideSettingsButton' src='/kernelspecs/polymake/close.svg' width=20px">
    		<img id='showSettingsButton_0' class='showSettingsButton' src='/kernelspecs/polymake/menu.svg' width=20px">
    <div id="model23894936808"></div>
    
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
       var container = document.getElementById( 'model23894936808' );
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
       var target = document.querySelector('#model23894936808');
    
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
    
       var objectnames = ["points","Vertices","unnamed__1","unnamed__2","unnamed__3","unnamed__4","unnamed__5"];
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 0, 0));
       allpoints.push(new PMPoint(2, 0, 0));
       allpoints.push(new PMPoint(1, 1, 0));
       allpoints.push(new PMPoint(0, 2, 0));
       allpoints.push(new PMPoint(2, 2, 0));
       allpoints.push(new PMPoint(1, 3, 0));
    
       <!-- Vertex style -->
       var points_material = new THREE.MeshBasicMaterial ( {color: 0x000000, } );
    
       points_material.side = THREE.DoubleSide;
       points_material.transparent = true;
    
       <!-- POINTS -->
       allpoints[0].makesphere(0.02,points_material);
       allpoints[1].makesphere(0.02,points_material);
       allpoints[2].makesphere(0.02,points_material);
       allpoints[3].makesphere(0.02,points_material);
       allpoints[4].makesphere(0.02,points_material);
       allpoints[5].makesphere(0.02,points_material);
       allpoints[0].makelabel("0");
       allpoints[1].makelabel("1");
       allpoints[2].makelabel("2");
       allpoints[3].makelabel("3");
       allpoints[4].makelabel("4");
       allpoints[5].makelabel("5");
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 0, 0));
       allpoints.push(new PMPoint(2, 0, 0));
       allpoints.push(new PMPoint(1, 1, 0));
       allpoints.push(new PMPoint(0, 2, 0));
       allpoints.push(new PMPoint(2, 2, 0));
       allpoints.push(new PMPoint(1, 3, 0));
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       <!-- Edge style -->
       var line_material = new THREE.LineBasicMaterial ( {color: 0x000000, linewidth: 1.5, } );
    
       line_material.side = THREE.DoubleSide;
       line_material.transparent = true;
    
       <!-- EDGES -->
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[3].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[3].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[4].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 0, 0));
       allpoints.push(new PMPoint(2, 0, 0));
       allpoints.push(new PMPoint(1, 1, 0));
    
       <!-- Vertex style -->
       var materials = [
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
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 0.5, side: THREE.DoubleSide , depthWrite: false, depthTest: false, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(2, 0, 1, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 0, 0));
       allpoints.push(new PMPoint(1, 1, 0));
       allpoints.push(new PMPoint(0, 2, 0));
    
       <!-- Vertex style -->
       var materials = [
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
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 0.5, side: THREE.DoubleSide , depthWrite: false, depthTest: false, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(2, 0, 1, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(1, 1, 0));
       allpoints.push(new PMPoint(0, 2, 0));
       allpoints.push(new PMPoint(2, 2, 0));
    
       <!-- Vertex style -->
       var materials = [
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
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 0.5, side: THREE.DoubleSide , depthWrite: false, depthTest: false, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(2, 1, 0, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 2, 0));
       allpoints.push(new PMPoint(2, 2, 0));
       allpoints.push(new PMPoint(1, 3, 0));
    
       <!-- Vertex style -->
       var materials = [
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
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 0.5, side: THREE.DoubleSide , depthWrite: false, depthTest: false, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(2, 0, 1, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(2, 0, 0));
       allpoints.push(new PMPoint(1, 1, 0));
       allpoints.push(new PMPoint(2, 2, 0));
    
       <!-- Vertex style -->
       var materials = [
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
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 0.5, side: THREE.DoubleSide , depthWrite: false, depthTest: false, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(2, 1, 0, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[2].vector);
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
    		document.getElementById('settings_0').style.position = 'absolute';
    		document.getElementById('settings_0').style.display = 'block';
    		document.getElementById('showSettingsButton_0').style.display = 'none';
    		document.getElementById('hideSettingsButton_0').style.display = 'block';
    		settingsShown = true;
    	}
    
    	function hideSettings(event){
    		event.currentTarget.style.display = 'none';
    		document.getElementById('settings_0').style.display = 'none';
    		document.getElementById('hideSettingsButton_0').style.display = 'none';
    		document.getElementById('showSettingsButton_0').style.display = 'block';
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
    
    		if (document.getElementById('tab_0').checked){
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
    							document.getElementById('transparency_0').style.display = 'block';
    							document.getElementById('transparencyRange_0').value = 1 - 
    								child.material.materials[k].opacity;
    							return;
    						}
    					}
    				} else if (	child.material.transparent && 
    								child.material.type == "MeshBasicMaterial" &&
    								child.geometry.type == "Geometry"){
    					document.getElementById('transparency_0').style.display = 'block';
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
    		document.getElementById('explodeRange_0').value = explodeScale;
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
    		document.getElementById('explodeRange_0').value = explodeScale;
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
    	var shownObjectsList = document.getElementById('shownObjectsList_0');
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
    tempobj = document.getElementById('explodeRange_0');
    if (tempobj) {
       tempobj.oninput = triggerExplode;
       document.getElementById('explodeCheckbox_0').onchange = triggerAutomaticExplode;
       document.getElementById('explodingSpeedRange_0').oninput = setExplodingSpeed;
    }
    tempobj = document.getElementById('foldRange_0');
    if (tempobj) {
       tempobj.oninput = fold;
    }
    document.getElementById('transparencyRange_0').oninput = changeTransparency;
    document.getElementById('changeRotationX_0').onchange = changeRotationX;
    document.getElementById('changeRotationY_0').onchange = changeRotationY;
    document.getElementById('changeRotationZ_0').onchange = changeRotationZ;
    document.getElementById('resetButton_0').onclick = resetScene;
    document.getElementById('rotationSpeedRange_0').oninput = changeRotationSpeedFactor;
    document.getElementById('labelsCheckboxInput_0').onchange = displayOrHideLabels;
    document.getElementById('takeScreenshot_0').onclick = takeSvgScreenshot;
    document.getElementById('showSettingsButton_0').onclick = showSettings;
    document.getElementById('hideSettingsButton_0').onclick = hideSettings;
    
    	
    
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
    		document.getElementById('hideSettingsButton_0').dispatchEvent(event);
    	} else{
    		document.getElementById('showSettingsButton_0').dispatchEvent(event);
    	}
    });
    
    if (foldable) moveToBaryCenter();
    
    
    });});
    // COMMON_CODE_BLOCK_END
    </script>
    
    </body>
    </html>



.. raw:: html

    <details><summary><pre style="display:inline"><small>Click here for additional output</small></pre></summary>
    <pre>
    polymake: used package threejs
       Three.js is a lightweight cross-browser JavaScript library/API used to create and display animated 3D computer graphics on a Web browser.
       See http://github.com/mrdoob for the source code.
    
    </pre>
    </details>




.. figure:: attachment:subdiv_of_points_3.gif
   :alt: {{:tutorial:subdiv_of_points_3.gif?300|}}

   {{:tutorial:subdiv_of_points_3.gif?300|}}

Indeed, the regular subdivision arising from this height function is a
triangulation. Since in this case point 2 is a vertex of the subdivision
it is colored red as well.

Generic Weights
~~~~~~~~~~~~~~~

The weights ``w`` are called *generic* if they define a triangulation of
the point set. If it is known in advance that the given weights are
generic, then we can obtain the subdivison as a simplicial complex. If
you want analyse your simplicial complex then you might want to use the
polymake objects: ``SimplicialComplex`` or
``GeometricSimplicialComplex``. The first one is a purely combinatoric
version. The latter one has an embedding into a space via the property
``COORDINATES``. Both objects live in the application ``topaz``. The
weights in the third example above are generic, so they define a
simplicial complex.


.. link

.. CODE-BLOCK:: perl

    polymake> $SC = new topaz::GeometricSimplicialComplex(COORDINATES=>$M, INPUT_FACES=>regular_subdivision($M,$w));

Note that ``regular_subdivision($M,$w)`` produces a list of the maximal
cells in the regular subdivision of the point configuration ``M``
respecting the weights ``w``. Now we can ask for many other properties,
e.g.


.. link

.. CODE-BLOCK:: perl

    polymake> print $SC->F_VECTOR;
    6 10 5




.. link

.. CODE-BLOCK:: perl

    polymake> print $SC->H_VECTOR;
    1 3 1 0




Check the properties of ``SimplicialComplex`` for more available
properties.

Non-generic Weights
~~~~~~~~~~~~~~~~~~~

If your weights are not generic or you are unsure, then you should
create a ``PolyhedralComplex``. Since every ``SubdivisionOfPoints``
Object has its underlying ``PolyhedralComplex`` as a property this can
be done simply by:


.. link

.. CODE-BLOCK:: perl

    polymake> $M = new Matrix<Rational>([[1,0,0,0],[1,0,1,0],[1,1,0,0],[1,1,1,0],[1,0,0,1],[1,1,0,1],[1,0,1,1],[1,1,1,1],[1,0,0,2]]);
    polymake> $w = new Vector<Rational>([1,0,0,1,0,1,1,0,1]);
    polymake> $S = new fan::SubdivisionOfPoints(POINTS=>$M,WEIGHTS=>$w);
    polymake> $PC = $S->POLYHEDRAL_COMPLEX;
    polymake> print $PC->MAXIMAL_POLYTOPES;
    {0 1 2 4}
    {1 2 4 7}
    {1 4 6 7 8}
    {2 4 5 7 8}
    {1 2 3 7}





::

   {0 1 2 4}
   {2 4 5 7 8}
   {1 2 3 7}
   {1 2 4 7}
   {1 4 6 7 8}

As with an simplicial complex you can do some computation with it. For
example:


.. link

.. CODE-BLOCK:: perl

    polymake> print $PC->F_VECTOR;
    9 20 17 5




.. link

.. CODE-BLOCK:: perl

    polymake> print $PC->SIMPLICIAL;
    false







.. link

.. CODE-BLOCK:: perl

    polymake> $PC->VISUAL;


.. raw:: html

    <!--
    polymake for andrew
    Thu Apr  4 10:32:51 2019
    pcom:
    -->
    
    
    <html>
       <head>
          <title>pcom:</title>
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
    
    		<div id='settings_1' class='settings'>
    			<div class=group id='explode_1'>
    				<strong>Explode</strong>
    				<input id='explodeRange_1' type='range' min=0 max=6 step=0.01 value=0>
    				<div class=indented><input id='explodeCheckbox_1' type='checkbox'>Automatic explosion</div>
    				<div class=suboption>Exploding speed</div>
    				<input id='explodingSpeedRange_1' type='range' min=0 max=0.5 step=0.001 value=0.05>
    			</div>
    
    			
    			<div class=group id='transparency_1' class='transparency'>
    				<strong>Transparency</strong>
    				<input id='transparencyRange_1' type='range' min=0 max=1 step=0.01 value=0>
    			</div>
    			
    			<div class=group id='rotation_1'>
    				<strong>Rotation</strong>
    				<div class=indented>
    					<div><input type='checkbox' id='changeRotationX_1'> x-axis</div>
    					<div><input type='checkbox' id='changeRotationY_1'> y-axis</div>
    					<div><input type='checkbox' id='changeRotationZ_1'> z-axis</div>
    					<button id='resetButton_1' class='resetButton' >Reset</button>
    				</div>
    
    				<div class=suboption>Rotation speed</div>
    				<input id='rotationSpeedRange_1' type='range' min=0 max=5 step=0.01 value=2>
    
    			</div>
    
    
    			<div class=group id='display_1'>
    				<strong>Display</strong>
    				<div class=indented>
    					<div id='shownObjectsList_1' class='shownObjectsList'></div>
    					<div class='labelsCheckbox'><input type='checkbox' id='labelsCheckboxInput_1' checked>Labels</div>
    				</div>
    			</div>
    
    
    			<div class=group id='svg_1'>
    				<strong>SVG</strong>
    				<div class=indented>
    					<form>
    						<input type="radio" name='screenshotMode' value='download' id='download_1' checked> Download<br>
    						<input type="radio" name='screenshotMode' value='tab' id='tab_1' > New tab<br>
    					</form>
    					<button id='takeScreenshot_1'>Screenshot</button>
    				</div>
    			</div>
    
    		</div>	<!-- end of settings -->
    		<img id='hideSettingsButton_1' style="display: none" class='hideSettingsButton' src='/kernelspecs/polymake/close.svg' width=20px">
    		<img id='showSettingsButton_1' class='showSettingsButton' src='/kernelspecs/polymake/menu.svg' width=20px">
    <div id="model24112417045"></div>
    
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
       var container = document.getElementById( 'model24112417045' );
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
       var target = document.querySelector('#model24112417045');
    
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
    
       var objectnames = ["unnamed__1","unnamed__2","unnamed__3","unnamed__4","unnamed__5"];
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 0, 0));
       allpoints.push(new PMPoint(0, 1, 0));
       allpoints.push(new PMPoint(1, 0, 0));
       allpoints.push(new PMPoint(0, 0, 1));
    
       <!-- Vertex style -->
       var points_material = new THREE.MeshBasicMaterial ( {color: 0xFF0000, } );
    
       points_material.side = THREE.DoubleSide;
       points_material.transparent = true;
    
       <!-- POINTS -->
       allpoints[0].makesphere(0.02,points_material);
       allpoints[1].makesphere(0.02,points_material);
       allpoints[2].makesphere(0.02,points_material);
       allpoints[3].makesphere(0.02,points_material);
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
       faces.vertices.push(allpoints[3].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(1, 3, 2, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(0, 3, 1, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(2, 3, 0, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(1, 2, 0, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 1, 0));
       allpoints.push(new PMPoint(1, 0, 0));
       allpoints.push(new PMPoint(0, 0, 1));
       allpoints.push(new PMPoint(1, 1, 1));
    
       <!-- Vertex style -->
       var points_material = new THREE.MeshBasicMaterial ( {color: 0xFF0000, } );
    
       points_material.side = THREE.DoubleSide;
       points_material.transparent = true;
    
       <!-- POINTS -->
       allpoints[0].makesphere(0.02,points_material);
       allpoints[1].makesphere(0.02,points_material);
       allpoints[2].makesphere(0.02,points_material);
       allpoints[3].makesphere(0.02,points_material);
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
       faces.vertices.push(allpoints[3].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(0, 3, 1, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(2, 3, 0, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(1, 3, 2, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(0, 1, 2, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 1, 0));
       allpoints.push(new PMPoint(0, 0, 1));
       allpoints.push(new PMPoint(0, 1, 1));
       allpoints.push(new PMPoint(1, 1, 1));
       allpoints.push(new PMPoint(0, 0, 2));
    
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
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(3, 1, 0, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(4, 1, 3, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(2, 4, 3, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(3, 0, 2, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(0, 1, 4, undefined, undefined, 0));
       faces.faces.push(new THREE.Face3(0, 4, 2, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[3].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[4].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[3].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(1, 0, 0));
       allpoints.push(new PMPoint(0, 0, 1));
       allpoints.push(new PMPoint(1, 0, 1));
       allpoints.push(new PMPoint(1, 1, 1));
       allpoints.push(new PMPoint(0, 0, 2));
    
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
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(1, 3, 0, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(1, 4, 3, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(2, 3, 4, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(0, 3, 2, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(1, 0, 2, undefined, undefined, 0));
       faces.faces.push(new THREE.Face3(1, 2, 4, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 1, 0));
       allpoints.push(new PMPoint(1, 0, 0));
       allpoints.push(new PMPoint(1, 1, 0));
       allpoints.push(new PMPoint(1, 1, 1));
    
       <!-- Vertex style -->
       var points_material = new THREE.MeshBasicMaterial ( {color: 0xFF0000, } );
    
       points_material.side = THREE.DoubleSide;
       points_material.transparent = true;
    
       <!-- POINTS -->
       allpoints[0].makesphere(0.02,points_material);
       allpoints[1].makesphere(0.02,points_material);
       allpoints[2].makesphere(0.02,points_material);
       allpoints[3].makesphere(0.02,points_material);
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
       faces.vertices.push(allpoints[3].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(2, 3, 1, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(0, 3, 2, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(1, 3, 0, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(2, 1, 0, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[2].vector);
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
    		document.getElementById('settings_1').style.position = 'absolute';
    		document.getElementById('settings_1').style.display = 'block';
    		document.getElementById('showSettingsButton_1').style.display = 'none';
    		document.getElementById('hideSettingsButton_1').style.display = 'block';
    		settingsShown = true;
    	}
    
    	function hideSettings(event){
    		event.currentTarget.style.display = 'none';
    		document.getElementById('settings_1').style.display = 'none';
    		document.getElementById('hideSettingsButton_1').style.display = 'none';
    		document.getElementById('showSettingsButton_1').style.display = 'block';
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
    
    		if (document.getElementById('tab_1').checked){
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
    							document.getElementById('transparency_1').style.display = 'block';
    							document.getElementById('transparencyRange_1').value = 1 - 
    								child.material.materials[k].opacity;
    							return;
    						}
    					}
    				} else if (	child.material.transparent && 
    								child.material.type == "MeshBasicMaterial" &&
    								child.geometry.type == "Geometry"){
    					document.getElementById('transparency_1').style.display = 'block';
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
    		document.getElementById('explodeRange_1').value = explodeScale;
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
    		document.getElementById('explodeRange_1').value = explodeScale;
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
    	var shownObjectsList = document.getElementById('shownObjectsList_1');
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
    tempobj = document.getElementById('explodeRange_1');
    if (tempobj) {
       tempobj.oninput = triggerExplode;
       document.getElementById('explodeCheckbox_1').onchange = triggerAutomaticExplode;
       document.getElementById('explodingSpeedRange_1').oninput = setExplodingSpeed;
    }
    tempobj = document.getElementById('foldRange_1');
    if (tempobj) {
       tempobj.oninput = fold;
    }
    document.getElementById('transparencyRange_1').oninput = changeTransparency;
    document.getElementById('changeRotationX_1').onchange = changeRotationX;
    document.getElementById('changeRotationY_1').onchange = changeRotationY;
    document.getElementById('changeRotationZ_1').onchange = changeRotationZ;
    document.getElementById('resetButton_1').onclick = resetScene;
    document.getElementById('rotationSpeedRange_1').oninput = changeRotationSpeedFactor;
    document.getElementById('labelsCheckboxInput_1').onchange = displayOrHideLabels;
    document.getElementById('takeScreenshot_1').onclick = takeSvgScreenshot;
    document.getElementById('showSettingsButton_1').onclick = showSettings;
    document.getElementById('hideSettingsButton_1').onclick = hideSettings;
    
    	
    
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
    		document.getElementById('hideSettingsButton_1').dispatchEvent(event);
    	} else{
    		document.getElementById('showSettingsButton_1').dispatchEvent(event);
    	}
    });
    
    if (foldable) moveToBaryCenter();
    
    
    });});
    // COMMON_CODE_BLOCK_END
    </script>
    
    </body>
    </html>



.. figure:: attachment:reg_subdiv_polytopal_subdiv.gif
   :alt: {{:tutorial:reg_subdiv_polytopal_subdiv.gif?300}}

   {{:tutorial:reg_subdiv_polytopal_subdiv.gif?300}}

Check the interactive help to get more information. For the image we
have used the javaview option *Explode Group of Geometries* to make the
cells of the subdivision visible.

Tropical Plücker Vectors and Matroid Decompositions of Hypersimplices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A *tropical Plücker vector* (which is a special lifting function on the
vertices of the *(d,n)*-hypersimplex induces a particularly interesting
kind of regular subdivision. The example below is for *d=2* and *n=4*.


.. link

.. CODE-BLOCK:: perl

    polymake> $p=new Vector<Rational>([1,0,0,0,0,0]);
    polymake> $msd=regular_subdivision(hypersimplex(2,4)->VERTICES,$p);

Each cell of this subdivision is a *matroid polytope*, that is, the
convex hull of characteristic vectors of the bases of some matroid on
*n* elements of rank *d*. The vertices of the hypersimplices happen to
be listed in lexicographical ordering. With this extra knowledge we can
cook up suitable labels for pretty-printing the maximal cells.


.. link

.. CODE-BLOCK:: perl

    polymake> print rows_labeled($msd,undef,["12","13","14","23","24","34"]);
    0:13 14 23 24 34
    1:12 13 14 23 24





In this case the *(2,4)*-hypersimplex (a.k.a. the regular octahedron) is
split into two egyptian pyramids.

Secondary cone and regularity of subdivisions
---------------------------------------------

The set of all the weight vectors which induce the same regular
subdivision form a convex cone. In this section we will learn how to
compute that cone with polymake and even how we can check whether a
given subdivision is regular or not.

Secondary Cone of a Regular Subdivision
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us look at the following regular subdivision (using homogeneous
coordinates):


.. link

.. CODE-BLOCK:: perl

    polymake> $points = new Matrix([[1,2,0],[1,0,2],[1,-2,0],[1,0,-2],[1,1,0],[1,0,1],[1,-1,0],[1,0,-1]]);
    polymake> $cells = new Array<Set<Int>>([[0,1,4,5],[0,3,4,7],[2,3,6,7],[1,2,5,6],[4,5,6,7]]);
    polymake> $S = new fan::SubdivisionOfPoints(POINTS=>$points, MAXIMAL_CELLS=>$cells);

we can visualize the subdivision with the following command:


.. link

.. CODE-BLOCK:: perl

    polymake> $S->VISUAL;


.. raw:: html

    <!--
    polymake for andrew
    Thu Apr  4 10:44:37 2019
    pcom:S
    -->
    
    
    <html>
       <head>
          <title>pcom:S</title>
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
    
    		<div id='settings_2' class='settings'>
    			<div class=group id='explode_2'>
    				<strong>Explode</strong>
    				<input id='explodeRange_2' type='range' min=0 max=6 step=0.01 value=0>
    				<div class=indented><input id='explodeCheckbox_2' type='checkbox'>Automatic explosion</div>
    				<div class=suboption>Exploding speed</div>
    				<input id='explodingSpeedRange_2' type='range' min=0 max=0.5 step=0.001 value=0.05>
    			</div>
    
    			
    			<div class=group id='transparency_2' class='transparency'>
    				<strong>Transparency</strong>
    				<input id='transparencyRange_2' type='range' min=0 max=1 step=0.01 value=0>
    			</div>
    			
    			<div class=group id='rotation_2'>
    				<strong>Rotation</strong>
    				<div class=indented>
    					<div><input type='checkbox' id='changeRotationX_2'> x-axis</div>
    					<div><input type='checkbox' id='changeRotationY_2'> y-axis</div>
    					<div><input type='checkbox' id='changeRotationZ_2'> z-axis</div>
    					<button id='resetButton_2' class='resetButton' >Reset</button>
    				</div>
    
    				<div class=suboption>Rotation speed</div>
    				<input id='rotationSpeedRange_2' type='range' min=0 max=5 step=0.01 value=2>
    
    			</div>
    
    
    			<div class=group id='display_2'>
    				<strong>Display</strong>
    				<div class=indented>
    					<div id='shownObjectsList_2' class='shownObjectsList'></div>
    					<div class='labelsCheckbox'><input type='checkbox' id='labelsCheckboxInput_2' checked>Labels</div>
    				</div>
    			</div>
    
    
    			<div class=group id='svg_2'>
    				<strong>SVG</strong>
    				<div class=indented>
    					<form>
    						<input type="radio" name='screenshotMode' value='download' id='download_2' checked> Download<br>
    						<input type="radio" name='screenshotMode' value='tab' id='tab_2' > New tab<br>
    					</form>
    					<button id='takeScreenshot_2'>Screenshot</button>
    				</div>
    			</div>
    
    		</div>	<!-- end of settings -->
    		<img id='hideSettingsButton_2' style="display: none" class='hideSettingsButton' src='/kernelspecs/polymake/close.svg' width=20px">
    		<img id='showSettingsButton_2' class='showSettingsButton' src='/kernelspecs/polymake/menu.svg' width=20px">
    <div id="model14313934235"></div>
    
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
       var container = document.getElementById( 'model14313934235' );
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
       var target = document.querySelector('#model14313934235');
    
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
    
       var objectnames = ["points","Vertices","unnamed__1","unnamed__2","unnamed__3","unnamed__4","unnamed__5"];
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(2, 0, 0));
       allpoints.push(new PMPoint(0, 2, 0));
       allpoints.push(new PMPoint(-2, 0, 0));
       allpoints.push(new PMPoint(0, -2, 0));
       allpoints.push(new PMPoint(1, 0, 0));
       allpoints.push(new PMPoint(0, 1, 0));
       allpoints.push(new PMPoint(-1, 0, 0));
       allpoints.push(new PMPoint(0, -1, 0));
    
       <!-- Vertex style -->
       var points_material = new THREE.MeshBasicMaterial ( {color: 0x000000, } );
    
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
       allpoints[7].makesphere(0.02,points_material);
       allpoints[0].makelabel("0");
       allpoints[1].makelabel("1");
       allpoints[2].makelabel("2");
       allpoints[3].makelabel("3");
       allpoints[4].makelabel("4");
       allpoints[5].makelabel("5");
       allpoints[6].makelabel("6");
       allpoints[7].makelabel("7");
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(2, 0, 0));
       allpoints.push(new PMPoint(0, 2, 0));
       allpoints.push(new PMPoint(-2, 0, 0));
       allpoints.push(new PMPoint(0, -2, 0));
       allpoints.push(new PMPoint(1, 0, 0));
       allpoints.push(new PMPoint(0, 1, 0));
       allpoints.push(new PMPoint(-1, 0, 0));
       allpoints.push(new PMPoint(0, -1, 0));
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       <!-- Edge style -->
       var line_material = new THREE.LineBasicMaterial ( {color: 0x000000, linewidth: 1.5, } );
    
       line_material.side = THREE.DoubleSide;
       line_material.transparent = true;
    
       <!-- EDGES -->
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[4].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[5].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[7].vector);
       line.vertices.push(allpoints[3].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[7].vector);
       line.vertices.push(allpoints[4].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[7].vector);
       line.vertices.push(allpoints[6].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(2, 0, 0));
       allpoints.push(new PMPoint(0, 2, 0));
       allpoints.push(new PMPoint(1, 0, 0));
       allpoints.push(new PMPoint(0, 1, 0));
    
       <!-- Vertex style -->
       var materials = [
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
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
       faces.vertices.push(allpoints[3].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 0.5, side: THREE.DoubleSide , depthWrite: false, depthTest: false, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(2, 0, 1, undefined, undefined, 0));
       faces.faces.push(new THREE.Face3(2, 1, 3, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(2, 0, 0));
       allpoints.push(new PMPoint(0, -2, 0));
       allpoints.push(new PMPoint(1, 0, 0));
       allpoints.push(new PMPoint(0, -1, 0));
    
       <!-- Vertex style -->
       var materials = [
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
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
       faces.vertices.push(allpoints[3].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 0.5, side: THREE.DoubleSide , depthWrite: false, depthTest: false, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(1, 0, 2, undefined, undefined, 0));
       faces.faces.push(new THREE.Face3(1, 2, 3, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(-2, 0, 0));
       allpoints.push(new PMPoint(0, -2, 0));
       allpoints.push(new PMPoint(-1, 0, 0));
       allpoints.push(new PMPoint(0, -1, 0));
    
       <!-- Vertex style -->
       var materials = [
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
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
       faces.vertices.push(allpoints[3].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 0.5, side: THREE.DoubleSide , depthWrite: false, depthTest: false, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(0, 1, 3, undefined, undefined, 0));
       faces.faces.push(new THREE.Face3(0, 3, 2, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 2, 0));
       allpoints.push(new PMPoint(-2, 0, 0));
       allpoints.push(new PMPoint(0, 1, 0));
       allpoints.push(new PMPoint(-1, 0, 0));
    
       <!-- Vertex style -->
       var materials = [
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
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
       faces.vertices.push(allpoints[3].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 0.5, side: THREE.DoubleSide , depthWrite: false, depthTest: false, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(3, 2, 0, undefined, undefined, 0));
       faces.faces.push(new THREE.Face3(3, 0, 1, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[3].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(1, 0, 0));
       allpoints.push(new PMPoint(0, 1, 0));
       allpoints.push(new PMPoint(-1, 0, 0));
       allpoints.push(new PMPoint(0, -1, 0));
    
       <!-- Vertex style -->
       var materials = [
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
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
       faces.vertices.push(allpoints[3].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 0.5, side: THREE.DoubleSide , depthWrite: false, depthTest: false, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(3, 0, 1, undefined, undefined, 0));
       faces.faces.push(new THREE.Face3(3, 1, 2, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[3].vector);
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
    		document.getElementById('settings_2').style.position = 'absolute';
    		document.getElementById('settings_2').style.display = 'block';
    		document.getElementById('showSettingsButton_2').style.display = 'none';
    		document.getElementById('hideSettingsButton_2').style.display = 'block';
    		settingsShown = true;
    	}
    
    	function hideSettings(event){
    		event.currentTarget.style.display = 'none';
    		document.getElementById('settings_2').style.display = 'none';
    		document.getElementById('hideSettingsButton_2').style.display = 'none';
    		document.getElementById('showSettingsButton_2').style.display = 'block';
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
    
    		if (document.getElementById('tab_2').checked){
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
    							document.getElementById('transparency_2').style.display = 'block';
    							document.getElementById('transparencyRange_2').value = 1 - 
    								child.material.materials[k].opacity;
    							return;
    						}
    					}
    				} else if (	child.material.transparent && 
    								child.material.type == "MeshBasicMaterial" &&
    								child.geometry.type == "Geometry"){
    					document.getElementById('transparency_2').style.display = 'block';
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
    		document.getElementById('explodeRange_2').value = explodeScale;
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
    		document.getElementById('explodeRange_2').value = explodeScale;
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
    	var shownObjectsList = document.getElementById('shownObjectsList_2');
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
    tempobj = document.getElementById('explodeRange_2');
    if (tempobj) {
       tempobj.oninput = triggerExplode;
       document.getElementById('explodeCheckbox_2').onchange = triggerAutomaticExplode;
       document.getElementById('explodingSpeedRange_2').oninput = setExplodingSpeed;
    }
    tempobj = document.getElementById('foldRange_2');
    if (tempobj) {
       tempobj.oninput = fold;
    }
    document.getElementById('transparencyRange_2').oninput = changeTransparency;
    document.getElementById('changeRotationX_2').onchange = changeRotationX;
    document.getElementById('changeRotationY_2').onchange = changeRotationY;
    document.getElementById('changeRotationZ_2').onchange = changeRotationZ;
    document.getElementById('resetButton_2').onclick = resetScene;
    document.getElementById('rotationSpeedRange_2').oninput = changeRotationSpeedFactor;
    document.getElementById('labelsCheckboxInput_2').onchange = displayOrHideLabels;
    document.getElementById('takeScreenshot_2').onclick = takeSvgScreenshot;
    document.getElementById('showSettingsButton_2').onclick = showSettings;
    document.getElementById('hideSettingsButton_2').onclick = hideSettings;
    
    	
    
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
    		document.getElementById('hideSettingsButton_2').dispatchEvent(event);
    	} else{
    		document.getElementById('showSettingsButton_2').dispatchEvent(event);
    	}
    });
    
    if (foldable) moveToBaryCenter();
    
    
    });});
    // COMMON_CODE_BLOCK_END
    </script>
    
    </body>
    </html>



.. figure:: attachment:square_reg_subdiv.png
   :alt: {{:tutorial:square_reg_subdiv.png?400|}}

   {{:tutorial:square_reg_subdiv.png?400|}}

As one can see from the picture this subdivision should be regular,
since we can easily find a weight vector which induces this subdivision.
Just lift all the points in the inner square to 0 and the points on the
outer square to 1. But now we want to take a look at all vectors which
induce this subdivision. This can be achieved by using the method
``secondary_cone``.


.. link

.. CODE-BLOCK:: perl

    polymake> $sc = $S->secondary_cone();
    polymake> print $sc->DIM;
    4




.. link

.. CODE-BLOCK:: perl

    polymake> print $sc->RAYS;
    1 1 1 1 0 0 0 0





.. link

.. CODE-BLOCK:: perl

    polymake> print $sc->LINEALITY_SPACE;
    3/2 3/2 -1/2 -1/2 1 1 0 0
    -13/14 15/14 23/14 -5/14 -2/7 5/7 1 0
    5/9 -1/3 5/9 13/9 5/9 1/9 5/9 1





We see that this cone is 4-dimensional with a 3-dimensional lineality
space. Every vector of the **interior** of that cone induces weights
that produces the same regular subdivision. The vectors on the boundary
give rise to a coarsening. As you see by looking at the coordinates, we
do not need to lift the inner square to 0, there are also weight vectors
which lift them in a different fashion but still produce the same
subdivision.

With the options ``lift_to_zero`` or ``lift_face_to_zero`` one can
specify points or a face which must be lifted to zero. Let us see what
happens if we force polymake to only allow weights that lift the inner
square to zero. As we see from our input, the face which describes the
inner square is the fifth one (but since polymake starts counting at
zero, it has number 4).


.. link

.. CODE-BLOCK:: perl

    polymake> $sc_fixed = $S->secondary_cone(lift_face_to_zero=>4);




.. link

.. CODE-BLOCK:: perl

    polymake> print $sc_fixed->DIM;
    1




.. link

.. CODE-BLOCK:: perl

    polymake> print $sc_fixed->RAYS;
    1 1 1 1 0 0 0 0





.. link

.. CODE-BLOCK:: perl

    polymake> print $sc_fixed->LINEALITY_SPACE;

Now you see that the lineality space of that cone is gone. The only way
to lift our points now is by raising all points on the outer square to
the same height which should be greater than 0.

Regularity of Subdivisions
~~~~~~~~~~~~~~~~~~~~~~~~~~

You might run into the situation where you want to check whether a
subdivision is regular or not. Don’t worry, polymake comes to rescue
along with it’s side-kick the ``is_regular`` function. The function
takes the same input as the ``secondary_cone`` function. The output is a
pair of a boolean and a vector. The boolean tells you whether or not the
subdivision is regular and the vector is a weight vector which induces
your subdivision. Let us check if our subdivision of the example above
is regular.


.. link

.. CODE-BLOCK:: perl

    polymake> print is_regular($points, $cells);
    1 <0 0 1 1 0 0 1/2 1/2>




.. link

.. CODE-BLOCK:: perl

    polymake> $lambda = new Vector("0 0 1 1 0 0 1/2 1/2");




.. link

.. CODE-BLOCK:: perl

    polymake> print regular_subdivision($points, $lambda);
    {4 5 6 7}
    {2 3 6 7}
    {1 2 5 6}
    {0 3 4 7}
    {0 1 4 5}





As you can see polymake tells us that our subdivision is indeed regular.
And the vector ``0 0 1 1 0 0 1/2 1/2`` produces the same subdivision we
started with. Now let us see what happens if we use a subdivision which
is not regular. We use the same points as above but different cells.


.. link

.. CODE-BLOCK:: perl

    polymake> $nreg_cells = new Array<Set<Int>>([[0,1,5],[0,4,5],[0,3,4],[3,4,7],[2,3,7],[2,6,7],[1,2,6],[1,5,6],[4,5,6,7]]);
    polymake> $nreg_pc = new fan::PolyhedralComplex(POINTS=>$points, INPUT_POLYTOPES=>$nreg_cells);
    polymake> $nreg_pc->VISUAL;


.. raw:: html

    <!--
    polymake for andrew
    Thu Apr  4 10:45:34 2019
    pcom:nreg_pc
    -->
    
    
    <html>
       <head>
          <title>pcom:nreg_pc</title>
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
    
    		<div id='settings_3' class='settings'>
    			<div class=group id='explode_3'>
    				<strong>Explode</strong>
    				<input id='explodeRange_3' type='range' min=0 max=6 step=0.01 value=0>
    				<div class=indented><input id='explodeCheckbox_3' type='checkbox'>Automatic explosion</div>
    				<div class=suboption>Exploding speed</div>
    				<input id='explodingSpeedRange_3' type='range' min=0 max=0.5 step=0.001 value=0.05>
    			</div>
    
    			
    			<div class=group id='transparency_3' class='transparency'>
    				<strong>Transparency</strong>
    				<input id='transparencyRange_3' type='range' min=0 max=1 step=0.01 value=0>
    			</div>
    			
    			<div class=group id='rotation_3'>
    				<strong>Rotation</strong>
    				<div class=indented>
    					<div><input type='checkbox' id='changeRotationX_3'> x-axis</div>
    					<div><input type='checkbox' id='changeRotationY_3'> y-axis</div>
    					<div><input type='checkbox' id='changeRotationZ_3'> z-axis</div>
    					<button id='resetButton_3' class='resetButton' >Reset</button>
    				</div>
    
    				<div class=suboption>Rotation speed</div>
    				<input id='rotationSpeedRange_3' type='range' min=0 max=5 step=0.01 value=2>
    
    			</div>
    
    
    			<div class=group id='display_3'>
    				<strong>Display</strong>
    				<div class=indented>
    					<div id='shownObjectsList_3' class='shownObjectsList'></div>
    					<div class='labelsCheckbox'><input type='checkbox' id='labelsCheckboxInput_3' checked>Labels</div>
    				</div>
    			</div>
    
    
    			<div class=group id='svg_3'>
    				<strong>SVG</strong>
    				<div class=indented>
    					<form>
    						<input type="radio" name='screenshotMode' value='download' id='download_3' checked> Download<br>
    						<input type="radio" name='screenshotMode' value='tab' id='tab_3' > New tab<br>
    					</form>
    					<button id='takeScreenshot_3'>Screenshot</button>
    				</div>
    			</div>
    
    		</div>	<!-- end of settings -->
    		<img id='hideSettingsButton_3' style="display: none" class='hideSettingsButton' src='/kernelspecs/polymake/close.svg' width=20px">
    		<img id='showSettingsButton_3' class='showSettingsButton' src='/kernelspecs/polymake/menu.svg' width=20px">
    <div id="model10016009662"></div>
    
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
       var container = document.getElementById( 'model10016009662' );
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
       var target = document.querySelector('#model10016009662');
    
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
    
       var objectnames = ["unnamed__1","unnamed__2","unnamed__3","unnamed__4","unnamed__5","unnamed__6","unnamed__7","unnamed__8","unnamed__9"];
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(2, 0, 0));
       allpoints.push(new PMPoint(0, 2, 0));
       allpoints.push(new PMPoint(0, 1, 0));
    
       <!-- Vertex style -->
       var points_material = new THREE.MeshBasicMaterial ( {color: 0xFF0000, } );
    
       points_material.side = THREE.DoubleSide;
       points_material.transparent = true;
    
       <!-- POINTS -->
       allpoints[0].makesphere(0.02,points_material);
       allpoints[1].makesphere(0.02,points_material);
       allpoints[2].makesphere(0.02,points_material);
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(1, 2, 0, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(2, 0, 0));
       allpoints.push(new PMPoint(0, 1, 0));
       allpoints.push(new PMPoint(1, 0, 0));
    
       <!-- Vertex style -->
       var points_material = new THREE.MeshBasicMaterial ( {color: 0xFF0000, } );
    
       points_material.side = THREE.DoubleSide;
       points_material.transparent = true;
    
       <!-- POINTS -->
       allpoints[0].makesphere(0.02,points_material);
       allpoints[1].makesphere(0.02,points_material);
       allpoints[2].makesphere(0.02,points_material);
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(0, 1, 2, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(2, 0, 0));
       allpoints.push(new PMPoint(1, 0, 0));
       allpoints.push(new PMPoint(0, -2, 0));
    
       <!-- Vertex style -->
       var points_material = new THREE.MeshBasicMaterial ( {color: 0xFF0000, } );
    
       points_material.side = THREE.DoubleSide;
       points_material.transparent = true;
    
       <!-- POINTS -->
       allpoints[0].makesphere(0.02,points_material);
       allpoints[1].makesphere(0.02,points_material);
       allpoints[2].makesphere(0.02,points_material);
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(1, 2, 0, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(1, 0, 0));
       allpoints.push(new PMPoint(0, -2, 0));
       allpoints.push(new PMPoint(0, -1, 0));
    
       <!-- Vertex style -->
       var points_material = new THREE.MeshBasicMaterial ( {color: 0xFF0000, } );
    
       points_material.side = THREE.DoubleSide;
       points_material.transparent = true;
    
       <!-- POINTS -->
       allpoints[0].makesphere(0.02,points_material);
       allpoints[1].makesphere(0.02,points_material);
       allpoints[2].makesphere(0.02,points_material);
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(2, 1, 0, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, -2, 0));
       allpoints.push(new PMPoint(0, -1, 0));
       allpoints.push(new PMPoint(-2, 0, 0));
    
       <!-- Vertex style -->
       var points_material = new THREE.MeshBasicMaterial ( {color: 0xFF0000, } );
    
       points_material.side = THREE.DoubleSide;
       points_material.transparent = true;
    
       <!-- POINTS -->
       allpoints[0].makesphere(0.02,points_material);
       allpoints[1].makesphere(0.02,points_material);
       allpoints[2].makesphere(0.02,points_material);
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(1, 2, 0, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, -1, 0));
       allpoints.push(new PMPoint(-2, 0, 0));
       allpoints.push(new PMPoint(-1, 0, 0));
    
       <!-- Vertex style -->
       var points_material = new THREE.MeshBasicMaterial ( {color: 0xFF0000, } );
    
       points_material.side = THREE.DoubleSide;
       points_material.transparent = true;
    
       <!-- POINTS -->
       allpoints[0].makesphere(0.02,points_material);
       allpoints[1].makesphere(0.02,points_material);
       allpoints[2].makesphere(0.02,points_material);
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(2, 1, 0, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 2, 0));
       allpoints.push(new PMPoint(-2, 0, 0));
       allpoints.push(new PMPoint(-1, 0, 0));
    
       <!-- Vertex style -->
       var points_material = new THREE.MeshBasicMaterial ( {color: 0xFF0000, } );
    
       points_material.side = THREE.DoubleSide;
       points_material.transparent = true;
    
       <!-- POINTS -->
       allpoints[0].makesphere(0.02,points_material);
       allpoints[1].makesphere(0.02,points_material);
       allpoints[2].makesphere(0.02,points_material);
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(0, 1, 2, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 2, 0));
       allpoints.push(new PMPoint(0, 1, 0));
       allpoints.push(new PMPoint(-1, 0, 0));
    
       <!-- Vertex style -->
       var points_material = new THREE.MeshBasicMaterial ( {color: 0xFF0000, } );
    
       points_material.side = THREE.DoubleSide;
       points_material.transparent = true;
    
       <!-- POINTS -->
       allpoints[0].makesphere(0.02,points_material);
       allpoints[1].makesphere(0.02,points_material);
       allpoints[2].makesphere(0.02,points_material);
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(0, 2, 1, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 1, 0));
       allpoints.push(new PMPoint(1, 0, 0));
       allpoints.push(new PMPoint(0, -1, 0));
       allpoints.push(new PMPoint(-1, 0, 0));
    
       <!-- Vertex style -->
       var points_material = new THREE.MeshBasicMaterial ( {color: 0xFF0000, } );
    
       points_material.side = THREE.DoubleSide;
       points_material.transparent = true;
    
       <!-- POINTS -->
       allpoints[0].makesphere(0.02,points_material);
       allpoints[1].makesphere(0.02,points_material);
       allpoints[2].makesphere(0.02,points_material);
       allpoints[3].makesphere(0.02,points_material);
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       var faces = new THREE.Geometry();
    
       <!-- VERTICES -->
       faces.vertices.push(allpoints[0].vector);
       faces.vertices.push(allpoints[1].vector);
       faces.vertices.push(allpoints[2].vector);
       faces.vertices.push(allpoints[3].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(2, 1, 0, undefined, undefined, 0));
       faces.faces.push(new THREE.Face3(2, 0, 3, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[2].vector);
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
    		document.getElementById('settings_3').style.position = 'absolute';
    		document.getElementById('settings_3').style.display = 'block';
    		document.getElementById('showSettingsButton_3').style.display = 'none';
    		document.getElementById('hideSettingsButton_3').style.display = 'block';
    		settingsShown = true;
    	}
    
    	function hideSettings(event){
    		event.currentTarget.style.display = 'none';
    		document.getElementById('settings_3').style.display = 'none';
    		document.getElementById('hideSettingsButton_3').style.display = 'none';
    		document.getElementById('showSettingsButton_3').style.display = 'block';
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
    
    		if (document.getElementById('tab_3').checked){
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
    							document.getElementById('transparency_3').style.display = 'block';
    							document.getElementById('transparencyRange_3').value = 1 - 
    								child.material.materials[k].opacity;
    							return;
    						}
    					}
    				} else if (	child.material.transparent && 
    								child.material.type == "MeshBasicMaterial" &&
    								child.geometry.type == "Geometry"){
    					document.getElementById('transparency_3').style.display = 'block';
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
    		document.getElementById('explodeRange_3').value = explodeScale;
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
    		document.getElementById('explodeRange_3').value = explodeScale;
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
    	var shownObjectsList = document.getElementById('shownObjectsList_3');
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
    tempobj = document.getElementById('explodeRange_3');
    if (tempobj) {
       tempobj.oninput = triggerExplode;
       document.getElementById('explodeCheckbox_3').onchange = triggerAutomaticExplode;
       document.getElementById('explodingSpeedRange_3').oninput = setExplodingSpeed;
    }
    tempobj = document.getElementById('foldRange_3');
    if (tempobj) {
       tempobj.oninput = fold;
    }
    document.getElementById('transparencyRange_3').oninput = changeTransparency;
    document.getElementById('changeRotationX_3').onchange = changeRotationX;
    document.getElementById('changeRotationY_3').onchange = changeRotationY;
    document.getElementById('changeRotationZ_3').onchange = changeRotationZ;
    document.getElementById('resetButton_3').onclick = resetScene;
    document.getElementById('rotationSpeedRange_3').oninput = changeRotationSpeedFactor;
    document.getElementById('labelsCheckboxInput_3').onchange = displayOrHideLabels;
    document.getElementById('takeScreenshot_3').onclick = takeSvgScreenshot;
    document.getElementById('showSettingsButton_3').onclick = showSettings;
    document.getElementById('hideSettingsButton_3').onclick = hideSettings;
    
    	
    
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
    		document.getElementById('hideSettingsButton_3').dispatchEvent(event);
    	} else{
    		document.getElementById('showSettingsButton_3').dispatchEvent(event);
    	}
    });
    
    if (foldable) moveToBaryCenter();
    
    
    });});
    // COMMON_CODE_BLOCK_END
    </script>
    
    </body>
    </html>



.. figure:: attachment:square_nreg_subdiv.png
   :alt: {{:tutorial:square_nreg_subdiv.png?400|}}

   {{:tutorial:square_nreg_subdiv.png?400|}}

This is quite similar to the ``mother of all examples`` (See:
“Triangulations” from DeLoera, Rambau and Santos). So with a similar
argument we can see that this subdivision is indeed not regular. How
does polymake feel about this?


.. link

.. CODE-BLOCK:: perl

    polymake> print is_regular($points, $nreg_cells);
    0 <>




Just as expected polymake tells us that the subdivision is not regular.
But although there is no weight function which induces this particular
subdivision we can still create a ``SubdivisionOfPoints`` object with
the given cells. What happens if we ask for the secondary cone of that
subdivision?


.. link

.. CODE-BLOCK:: perl

    polymake> $nreg_sub = new fan::SubdivisionOfPoints(POINTS=>$points, MAXIMAL_CELLS=>$nreg_cells);
    polymake> $nreg_sc = $nreg_sub->secondary_cone();
    polymake> print $nreg_sc->DIM;
    4




This might be weird at first. Although our subdivision is not regular
(meaning there does not exist a secondary cone for that subdivision) the
secondary cone seems to be 4 dimensional. Let us investigate.


.. link

.. CODE-BLOCK:: perl

    polymake> print $nreg_sc->RAYS;
    1 1 1 1 0 0 0 0





.. link

.. CODE-BLOCK:: perl

    polymake> print $nreg_sc->LINEALITY_SPACE;
    3/2 3/2 -1/2 -1/2 1 1 0 0
    -13/14 15/14 23/14 -5/14 -2/7 5/7 1 0
    5/9 -1/3 5/9 13/9 5/9 1/9 5/9 1





This looks familiar, doesn’t it? We got the same secondary cone as we
did in our first example. The reason for this is that the first example
is the finest regular coarsening of our non regular subdivision. And
this is exactly what we see here. In other words, the subdivision of our
first example is the closest regular thing similar to our non regular
subdivision and this is exactly what the ``secondary_cone`` function
produces. If a subdivision is not regular it gives you the secondary
cone of the next best regular subdivision.
