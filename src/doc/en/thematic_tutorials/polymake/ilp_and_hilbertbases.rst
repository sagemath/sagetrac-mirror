.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


ILP and Hilbert bases
---------------------

A first example
~~~~~~~~~~~~~~~

First we will construct a new rational polytope:


::

    polymake> $p=new Polytope<Rational>;




::

    polymake> $p->POINTS=<<".";
    ........> 1 0 0 0
    ........> 1 1 0 0
    ........> 1 0 1 0
    ........> 1 1 1 0
    ........> 1 0 0 1
    ........> 1 1 0 1
    ........> 1 0 1 1
    ........> 1 1 1 1
    ........> .

Note that points in ``polymake`` are always given in homogenous
coordinates. I.e., the point (a,b,c) in R3 is represented as ``1 a b c``
in ``polymake``.

Now we can examine some properties of ``$p``. For instance we can
determine the number of facets or whether ``$p`` is simple:


::

    polymake> print $p->N_FACETS;
    6




.. raw:: html

    <details><summary><pre style="display:inline"><small>Click here for additional output</small></pre></summary>
    <pre>
    polymake: used package cdd
      cddlib
      Implementation of the double description method of Motzkin et al.
      Copyright by Komei Fukuda.
      http://www-oldurls.inf.ethz.ch/personal/fukudak/cdd_home/
    
    polymake: used package lrs
      Implementation of the reverse search algorithm of Avis and Fukuda.
      Copyright by David Avis.
      http://cgm.cs.mcgill.ca/~avis/C/lrs.html
    
    </pre>
    </details>




::

    polymake> print $p->SIMPLE;
    true




As you might already have noticed, our polytope is just a 3-dimensional
cube. So there would have been an easier way to create it using the
client ``cube``:


::

    polymake> $c = cube(3,0);

(You can check out the details of any function in the ```polymake``
documentation <http://wwwopt.mathematik.tu-darmstadt.de/polymake_doku/2.9.8/>`__.)

And we can also verify that the two polytopes are actually equal:


::

    polymake> print equal_polyhedra($p,$c);
    true




Another example
~~~~~~~~~~~~~~~

Now let us proceed with a somewhat more interesting example: The convex
hull of 20 randomly chosen points on the 2-dimensional sphere. |{{
:tutorial:ilp:rand_sphere.png?200|}}|

.. |{{ :tutorial:ilp:rand_sphere.png?200|}}| image:: attachment:rand_sphere.png


::

    polymake> $rs = rand_sphere(3,20);

``polymake`` can of course visualise this polytope:


::

    polymake> $rs->VISUAL;


.. raw:: html

    <!--
    polymake for andrew
    Thu Mar 28 11:23:22 2019
    rs
    -->
    
    
    <html>
       <head>
          <title>rs</title>
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
    <div id="model38065534069"></div>
    
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
       var container = document.getElementById( 'model38065534069' );
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
       var target = document.querySelector('#model38065534069');
    
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
    
       var objectnames = ["rs"];
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0.0347759, -0.517072, -0.855235));
       allpoints.push(new PMPoint(0.0940138, -0.463518, -0.881086));
       allpoints.push(new PMPoint(-0.0102355, -0.492099, 0.870479));
       allpoints.push(new PMPoint(-0.303517, 0.865815, 0.397797));
       allpoints.push(new PMPoint(0.660284, 0.202359, -0.72324));
       allpoints.push(new PMPoint(-0.102607, -0.968554, 0.226659));
       allpoints.push(new PMPoint(-0.993138, -0.090402, -0.0741867));
       allpoints.push(new PMPoint(-0.660087, 0.374854, 0.650976));
       allpoints.push(new PMPoint(-0.0652556, 0.916093, -0.395619));
       allpoints.push(new PMPoint(-0.673716, 0.464125, -0.575061));
       allpoints.push(new PMPoint(-0.304672, -0.776212, 0.55197));
       allpoints.push(new PMPoint(0.304639, 0.436836, -0.846386));
       allpoints.push(new PMPoint(0.736753, 0.500536, 0.454597));
       allpoints.push(new PMPoint(-0.703485, 0.679363, -0.208746));
       allpoints.push(new PMPoint(0.72315, -0.163198, -0.671133));
       allpoints.push(new PMPoint(0.98943, -0.140866, -0.034413));
       allpoints.push(new PMPoint(0.742389, 0.440088, 0.505155));
       allpoints.push(new PMPoint(-0.0898398, 0.96969, -0.227221));
       allpoints.push(new PMPoint(-0.00629402, -0.178753, 0.983874));
       allpoints.push(new PMPoint(-0.382233, -0.568978, -0.728123));
    
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
       allpoints[7].makesphere(0.02,points_material);
       allpoints[8].makesphere(0.02,points_material);
       allpoints[9].makesphere(0.02,points_material);
       allpoints[10].makesphere(0.02,points_material);
       allpoints[11].makesphere(0.02,points_material);
       allpoints[12].makesphere(0.02,points_material);
       allpoints[13].makesphere(0.02,points_material);
       allpoints[14].makesphere(0.02,points_material);
       allpoints[15].makesphere(0.02,points_material);
       allpoints[16].makesphere(0.02,points_material);
       allpoints[17].makesphere(0.02,points_material);
       allpoints[18].makesphere(0.02,points_material);
       allpoints[19].makesphere(0.02,points_material);
       allpoints[0].makelabel("0");
       allpoints[1].makelabel("1");
       allpoints[2].makelabel("2");
       allpoints[3].makelabel("3");
       allpoints[4].makelabel("4");
       allpoints[5].makelabel("5");
       allpoints[6].makelabel("6");
       allpoints[7].makelabel("7");
       allpoints[8].makelabel("8");
       allpoints[9].makelabel("9");
       allpoints[10].makelabel("10");
       allpoints[11].makelabel("11");
       allpoints[12].makelabel("12");
       allpoints[13].makelabel("13");
       allpoints[14].makelabel("14");
       allpoints[15].makelabel("15");
       allpoints[16].makelabel("16");
       allpoints[17].makelabel("17");
       allpoints[18].makelabel("18");
       allpoints[19].makelabel("19");
    
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
       faces.vertices.push(allpoints[7].vector);
       faces.vertices.push(allpoints[8].vector);
       faces.vertices.push(allpoints[9].vector);
       faces.vertices.push(allpoints[10].vector);
       faces.vertices.push(allpoints[11].vector);
       faces.vertices.push(allpoints[12].vector);
       faces.vertices.push(allpoints[13].vector);
       faces.vertices.push(allpoints[14].vector);
       faces.vertices.push(allpoints[15].vector);
       faces.vertices.push(allpoints[16].vector);
       faces.vertices.push(allpoints[17].vector);
       faces.vertices.push(allpoints[18].vector);
       faces.vertices.push(allpoints[19].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var faces_material = new THREE.MeshBasicMaterial ( {color: 0x77EC9E, transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, } );
    
       faces_material.side = THREE.DoubleSide;
       faces_material.transparent = true;
    
       <!-- FACETS --> 
       faces.faces.push(new THREE.Face3(3, 16, 12, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(18, 16, 3, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(2, 16, 18, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(4, 14, 1, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(1, 14, 0, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(14, 5, 0, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(11, 1, 9, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(8, 11, 9, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(0, 5, 19, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(5, 2, 10, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(3, 17, 13, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(3, 13, 7, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(6, 9, 19, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(7, 13, 6, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(13, 9, 6, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(7, 6, 10, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(5, 6, 19, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(5, 10, 6, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(13, 8, 9, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(17, 8, 13, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(3, 7, 18, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(10, 18, 7, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(2, 18, 10, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(1, 19, 9, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(1, 0, 19, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(8, 4, 11, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(4, 1, 11, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(17, 12, 8, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(3, 12, 17, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(12, 4, 8, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(4, 15, 14, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(12, 15, 4, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(15, 5, 14, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(15, 2, 5, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(15, 16, 2, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(12, 16, 15, undefined, undefined, 0));
    
    
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
       line.vertices.push(allpoints[16].vector);
       line.vertices.push(allpoints[12].vector);
       line.vertices.push(allpoints[3].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[18].vector);
       line.vertices.push(allpoints[16].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[18].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[16].vector);
       line.vertices.push(allpoints[18].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[14].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[4].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[14].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[14].vector);
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[14].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[11].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[9].vector);
       line.vertices.push(allpoints[11].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[8].vector);
       line.vertices.push(allpoints[11].vector);
       line.vertices.push(allpoints[9].vector);
       line.vertices.push(allpoints[8].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[19].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[10].vector);
       line.vertices.push(allpoints[5].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[17].vector);
       line.vertices.push(allpoints[13].vector);
       line.vertices.push(allpoints[3].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[13].vector);
       line.vertices.push(allpoints[7].vector);
       line.vertices.push(allpoints[3].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[9].vector);
       line.vertices.push(allpoints[19].vector);
       line.vertices.push(allpoints[6].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[7].vector);
       line.vertices.push(allpoints[13].vector);
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[7].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[13].vector);
       line.vertices.push(allpoints[9].vector);
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[13].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[7].vector);
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[10].vector);
       line.vertices.push(allpoints[7].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[19].vector);
       line.vertices.push(allpoints[5].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[10].vector);
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[5].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[13].vector);
       line.vertices.push(allpoints[8].vector);
       line.vertices.push(allpoints[9].vector);
       line.vertices.push(allpoints[13].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[17].vector);
       line.vertices.push(allpoints[8].vector);
       line.vertices.push(allpoints[13].vector);
       line.vertices.push(allpoints[17].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[7].vector);
       line.vertices.push(allpoints[18].vector);
       line.vertices.push(allpoints[3].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[10].vector);
       line.vertices.push(allpoints[18].vector);
       line.vertices.push(allpoints[7].vector);
       line.vertices.push(allpoints[10].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[18].vector);
       line.vertices.push(allpoints[10].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[19].vector);
       line.vertices.push(allpoints[9].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[19].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[8].vector);
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[11].vector);
       line.vertices.push(allpoints[8].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[11].vector);
       line.vertices.push(allpoints[4].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[17].vector);
       line.vertices.push(allpoints[12].vector);
       line.vertices.push(allpoints[8].vector);
       line.vertices.push(allpoints[17].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[12].vector);
       line.vertices.push(allpoints[17].vector);
       line.vertices.push(allpoints[3].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[12].vector);
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[8].vector);
       line.vertices.push(allpoints[12].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[15].vector);
       line.vertices.push(allpoints[14].vector);
       line.vertices.push(allpoints[4].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[12].vector);
       line.vertices.push(allpoints[15].vector);
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[12].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[15].vector);
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[14].vector);
       line.vertices.push(allpoints[15].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[15].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[15].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[15].vector);
       line.vertices.push(allpoints[16].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[15].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[12].vector);
       line.vertices.push(allpoints[16].vector);
       line.vertices.push(allpoints[15].vector);
       line.vertices.push(allpoints[12].vector);
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




Now we will create yet another new polytope by scaling our random sphere
by a factor lambda. (Otherwise there are rather few integral points
contained in it.)

To this end, we have to multiply every coordinate (except for the
homogenising 1 in the beginning) of every vertex by lamda. Then we can
create a new polytope by specifying its vertices.

.. figure:: attachment:rand_sphere_lattice.png
   :alt: {{ :tutorial:ilp:rand_sphere_lattice.png?200|}}

   {{ :tutorial:ilp:rand_sphere_lattice.png?200|}}


::

    polymake> $lambda=2;




::

    polymake> $s=new Matrix<Rational>([[1,0,0,0],[0,$lambda,0,0],[0,0,$lambda,0],[0,0,0,$lambda]]);




::

    polymake> print $s;
    1 0 0 0
    0 2 0 0
    0 0 2 0
    0 0 0 2





::

    polymake> $scaled_rs=new Polytope<Rational>(VERTICES=>($rs->VERTICES * $s), LINEALITY_SPACE=>[]);

``polymake`` can visualise the polytope together with its lattice
points:


::

    polymake> $scaled_rs->VISUAL->LATTICE_COLORED;


.. raw:: html

    <!--
    polymake for andrew
    Thu Mar 28 11:23:37 2019
    scaled_rs
    -->
    
    
    <html>
       <head>
          <title>scaled_rs</title>
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
    <div id="model2138786213"></div>
    
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
       var container = document.getElementById( 'model2138786213' );
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
       var target = document.querySelector('#model2138786213');
    
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
    
       var objectnames = ["scaled_rs","Lattice points and vertices of scaled_rs"];
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0.0695518, -1.03414, -1.71047));
       allpoints.push(new PMPoint(0.188028, -0.927036, -1.76217));
       allpoints.push(new PMPoint(-0.020471, -0.984198, 1.74096));
       allpoints.push(new PMPoint(-0.607035, 1.73163, 0.795593));
       allpoints.push(new PMPoint(1.32057, 0.404718, -1.44648));
       allpoints.push(new PMPoint(-0.205214, -1.93711, 0.453317));
       allpoints.push(new PMPoint(-1.98628, -0.180804, -0.148373));
       allpoints.push(new PMPoint(-1.32017, 0.749708, 1.30195));
       allpoints.push(new PMPoint(-0.130511, 1.83219, -0.791238));
       allpoints.push(new PMPoint(-1.34743, 0.928249, -1.15012));
       allpoints.push(new PMPoint(-0.609343, -1.55242, 1.10394));
       allpoints.push(new PMPoint(0.609278, 0.873673, -1.69277));
       allpoints.push(new PMPoint(1.47351, 1.00107, 0.909194));
       allpoints.push(new PMPoint(-1.40697, 1.35873, -0.417491));
       allpoints.push(new PMPoint(1.4463, -0.326397, -1.34227));
       allpoints.push(new PMPoint(1.97886, -0.281731, -0.0688259));
       allpoints.push(new PMPoint(1.48478, 0.880175, 1.01031));
       allpoints.push(new PMPoint(-0.17968, 1.93938, -0.454442));
       allpoints.push(new PMPoint(-0.012588, -0.357506, 1.96775));
       allpoints.push(new PMPoint(-0.764466, -1.13796, -1.45625));
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       <!-- Edge style -->
       var line_material = new THREE.LineBasicMaterial ( {color: 0x000000, linewidth: 1.5, } );
    
       line_material.side = THREE.DoubleSide;
       line_material.transparent = true;
    
       <!-- EDGES --> 
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[16].vector);
       line.vertices.push(allpoints[12].vector);
       line.vertices.push(allpoints[3].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[18].vector);
       line.vertices.push(allpoints[16].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[18].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[16].vector);
       line.vertices.push(allpoints[18].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[14].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[4].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[14].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[14].vector);
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[14].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[11].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[9].vector);
       line.vertices.push(allpoints[11].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[8].vector);
       line.vertices.push(allpoints[11].vector);
       line.vertices.push(allpoints[9].vector);
       line.vertices.push(allpoints[8].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[19].vector);
       line.vertices.push(allpoints[0].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[10].vector);
       line.vertices.push(allpoints[5].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[17].vector);
       line.vertices.push(allpoints[13].vector);
       line.vertices.push(allpoints[3].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[13].vector);
       line.vertices.push(allpoints[7].vector);
       line.vertices.push(allpoints[3].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[9].vector);
       line.vertices.push(allpoints[19].vector);
       line.vertices.push(allpoints[6].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[7].vector);
       line.vertices.push(allpoints[13].vector);
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[7].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[13].vector);
       line.vertices.push(allpoints[9].vector);
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[13].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[7].vector);
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[10].vector);
       line.vertices.push(allpoints[7].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[19].vector);
       line.vertices.push(allpoints[5].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[10].vector);
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[5].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[13].vector);
       line.vertices.push(allpoints[8].vector);
       line.vertices.push(allpoints[9].vector);
       line.vertices.push(allpoints[13].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[17].vector);
       line.vertices.push(allpoints[8].vector);
       line.vertices.push(allpoints[13].vector);
       line.vertices.push(allpoints[17].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[7].vector);
       line.vertices.push(allpoints[18].vector);
       line.vertices.push(allpoints[3].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[10].vector);
       line.vertices.push(allpoints[18].vector);
       line.vertices.push(allpoints[7].vector);
       line.vertices.push(allpoints[10].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[18].vector);
       line.vertices.push(allpoints[10].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[19].vector);
       line.vertices.push(allpoints[9].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[19].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[8].vector);
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[11].vector);
       line.vertices.push(allpoints[8].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[11].vector);
       line.vertices.push(allpoints[4].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[17].vector);
       line.vertices.push(allpoints[12].vector);
       line.vertices.push(allpoints[8].vector);
       line.vertices.push(allpoints[17].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[12].vector);
       line.vertices.push(allpoints[17].vector);
       line.vertices.push(allpoints[3].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[12].vector);
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[8].vector);
       line.vertices.push(allpoints[12].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[15].vector);
       line.vertices.push(allpoints[14].vector);
       line.vertices.push(allpoints[4].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[12].vector);
       line.vertices.push(allpoints[15].vector);
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[12].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[15].vector);
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[14].vector);
       line.vertices.push(allpoints[15].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[15].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[15].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[15].vector);
       line.vertices.push(allpoints[16].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[15].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[12].vector);
       line.vertices.push(allpoints[16].vector);
       line.vertices.push(allpoints[15].vector);
       line.vertices.push(allpoints[12].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(-1, -1, 0));
       allpoints.push(new PMPoint(-1, 0, -1));
       allpoints.push(new PMPoint(-1, 0, 0));
       allpoints.push(new PMPoint(-1, 0, 1));
       allpoints.push(new PMPoint(-1, 1, -1));
       allpoints.push(new PMPoint(-1, 1, 0));
       allpoints.push(new PMPoint(-1, 1, 1));
       allpoints.push(new PMPoint(0, -1, -1));
       allpoints.push(new PMPoint(0, -1, 0));
       allpoints.push(new PMPoint(0, -1, 1));
       allpoints.push(new PMPoint(0, 0, -1));
       allpoints.push(new PMPoint(0, 0, 0));
       allpoints.push(new PMPoint(0, 0, 1));
       allpoints.push(new PMPoint(0, 1, -1));
       allpoints.push(new PMPoint(0, 1, 0));
       allpoints.push(new PMPoint(0, 1, 1));
       allpoints.push(new PMPoint(1, 0, -1));
       allpoints.push(new PMPoint(1, 0, 0));
       allpoints.push(new PMPoint(1, 0, 1));
       allpoints.push(new PMPoint(1, 1, 0));
       allpoints.push(new PMPoint(1, 1, 1));
       allpoints.push(new PMPoint(0.0695518, -1.03414, -1.71047));
       allpoints.push(new PMPoint(0.188028, -0.927036, -1.76217));
       allpoints.push(new PMPoint(-0.020471, -0.984198, 1.74096));
       allpoints.push(new PMPoint(-0.607035, 1.73163, 0.795593));
       allpoints.push(new PMPoint(1.32057, 0.404718, -1.44648));
       allpoints.push(new PMPoint(-0.205214, -1.93711, 0.453317));
       allpoints.push(new PMPoint(-1.98628, -0.180804, -0.148373));
       allpoints.push(new PMPoint(-1.32017, 0.749708, 1.30195));
       allpoints.push(new PMPoint(-0.130511, 1.83219, -0.791238));
       allpoints.push(new PMPoint(-1.34743, 0.928249, -1.15012));
       allpoints.push(new PMPoint(-0.609343, -1.55242, 1.10394));
       allpoints.push(new PMPoint(0.609278, 0.873673, -1.69277));
       allpoints.push(new PMPoint(1.47351, 1.00107, 0.909194));
       allpoints.push(new PMPoint(-1.40697, 1.35873, -0.417491));
       allpoints.push(new PMPoint(1.4463, -0.326397, -1.34227));
       allpoints.push(new PMPoint(1.97886, -0.281731, -0.0688259));
       allpoints.push(new PMPoint(1.48478, 0.880175, 1.01031));
       allpoints.push(new PMPoint(-0.17968, 1.93938, -0.454442));
       allpoints.push(new PMPoint(-0.012588, -0.357506, 1.96775));
       allpoints.push(new PMPoint(-0.764466, -1.13796, -1.45625));
    
       <!-- Vertex style -->
       var materials = [
          new THREE.MeshBasicMaterial({ color: 0x1EFA1E, }),
          new THREE.MeshBasicMaterial({ color: 0x1EFA1E, }),
          new THREE.MeshBasicMaterial({ color: 0x1EFA1E, }),
          new THREE.MeshBasicMaterial({ color: 0x1EFA1E, }),
          new THREE.MeshBasicMaterial({ color: 0x1EFA1E, }),
          new THREE.MeshBasicMaterial({ color: 0x1EFA1E, }),
          new THREE.MeshBasicMaterial({ color: 0x1EFA1E, }),
          new THREE.MeshBasicMaterial({ color: 0x1EFA1E, }),
          new THREE.MeshBasicMaterial({ color: 0x1EFA1E, }),
          new THREE.MeshBasicMaterial({ color: 0x1EFA1E, }),
          new THREE.MeshBasicMaterial({ color: 0x1EFA1E, }),
          new THREE.MeshBasicMaterial({ color: 0x1EFA1E, }),
          new THREE.MeshBasicMaterial({ color: 0x1EFA1E, }),
          new THREE.MeshBasicMaterial({ color: 0x1EFA1E, }),
          new THREE.MeshBasicMaterial({ color: 0x1EFA1E, }),
          new THREE.MeshBasicMaterial({ color: 0x1EFA1E, }),
          new THREE.MeshBasicMaterial({ color: 0x1EFA1E, }),
          new THREE.MeshBasicMaterial({ color: 0x1EFA1E, }),
          new THREE.MeshBasicMaterial({ color: 0x1EFA1E, }),
          new THREE.MeshBasicMaterial({ color: 0x1EFA1E, }),
          new THREE.MeshBasicMaterial({ color: 0x1EFA1E, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
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
       allpoints[7].makesphere(0.02,materials[7]);
       allpoints[8].makesphere(0.02,materials[8]);
       allpoints[9].makesphere(0.02,materials[9]);
       allpoints[10].makesphere(0.02,materials[10]);
       allpoints[11].makesphere(0.02,materials[11]);
       allpoints[12].makesphere(0.02,materials[12]);
       allpoints[13].makesphere(0.02,materials[13]);
       allpoints[14].makesphere(0.02,materials[14]);
       allpoints[15].makesphere(0.02,materials[15]);
       allpoints[16].makesphere(0.02,materials[16]);
       allpoints[17].makesphere(0.02,materials[17]);
       allpoints[18].makesphere(0.02,materials[18]);
       allpoints[19].makesphere(0.02,materials[19]);
       allpoints[20].makesphere(0.02,materials[20]);
       allpoints[21].makesphere(0.02,materials[21]);
       allpoints[22].makesphere(0.02,materials[22]);
       allpoints[23].makesphere(0.02,materials[23]);
       allpoints[24].makesphere(0.02,materials[24]);
       allpoints[25].makesphere(0.02,materials[25]);
       allpoints[26].makesphere(0.02,materials[26]);
       allpoints[27].makesphere(0.02,materials[27]);
       allpoints[28].makesphere(0.02,materials[28]);
       allpoints[29].makesphere(0.02,materials[29]);
       allpoints[30].makesphere(0.02,materials[30]);
       allpoints[31].makesphere(0.02,materials[31]);
       allpoints[32].makesphere(0.02,materials[32]);
       allpoints[33].makesphere(0.02,materials[33]);
       allpoints[34].makesphere(0.02,materials[34]);
       allpoints[35].makesphere(0.02,materials[35]);
       allpoints[36].makesphere(0.02,materials[36]);
       allpoints[37].makesphere(0.02,materials[37]);
       allpoints[38].makesphere(0.02,materials[38]);
       allpoints[39].makesphere(0.02,materials[39]);
       allpoints[40].makesphere(0.02,materials[40]);
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
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



Now will construct the integer hull of ``$scaled_rs`` and visualise it:

.. figure:: attachment:ilp_lattice.png
   :alt: {{ :tutorial:ilp:ilp_lattice.png?200|}}

   {{ :tutorial:ilp:ilp_lattice.png?200|}}


::

    polymake> $integer_hull=new Polytope<Rational>(POINTS=>$scaled_rs->LATTICE_POINTS);




::

    polymake> $integer_hull->VISUAL->LATTICE_COLORED;


.. raw:: html

    <!--
    polymake for andrew
    Thu Mar 28 11:23:42 2019
    integer_hull
    -->
    
    
    <html>
       <head>
          <title>integer_hull</title>
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
    <div id="model64990790564"></div>
    
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
       var container = document.getElementById( 'model64990790564' );
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
       var target = document.querySelector('#model64990790564');
    
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
    
       var objectnames = ["integer_hull","Lattice points and vertices of integer_hull"];
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(-1, -1, 0));
       allpoints.push(new PMPoint(-1, 0, -1));
       allpoints.push(new PMPoint(-1, 0, 1));
       allpoints.push(new PMPoint(-1, 1, -1));
       allpoints.push(new PMPoint(-1, 1, 1));
       allpoints.push(new PMPoint(0, -1, -1));
       allpoints.push(new PMPoint(0, -1, 1));
       allpoints.push(new PMPoint(0, 1, -1));
       allpoints.push(new PMPoint(1, 0, -1));
       allpoints.push(new PMPoint(1, 0, 1));
       allpoints.push(new PMPoint(1, 1, 0));
       allpoints.push(new PMPoint(1, 1, 1));
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
       <!-- Edge style -->
       var line_material = new THREE.LineBasicMaterial ( {color: 0x000000, linewidth: 1.5, } );
    
       line_material.side = THREE.DoubleSide;
       line_material.transparent = true;
    
       <!-- EDGES --> 
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[8].vector);
       line.vertices.push(allpoints[7].vector);
       line.vertices.push(allpoints[10].vector);
       line.vertices.push(allpoints[8].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[7].vector);
       line.vertices.push(allpoints[8].vector);
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[3].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[5].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[4].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[8].vector);
       line.vertices.push(allpoints[9].vector);
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[8].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[11].vector);
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[9].vector);
       line.vertices.push(allpoints[11].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[8].vector);
       line.vertices.push(allpoints[10].vector);
       line.vertices.push(allpoints[11].vector);
       line.vertices.push(allpoints[9].vector);
       line.vertices.push(allpoints[8].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[10].vector);
       line.vertices.push(allpoints[7].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[11].vector);
       line.vertices.push(allpoints[10].vector);
       obj.add(new THREE.Line(line, line_material));
    
       scene.add(obj);
       all_objects.push(obj);
    
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(0, 0, 0));
       allpoints.push(new PMPoint(-1, -1, 0));
       allpoints.push(new PMPoint(-1, 0, -1));
       allpoints.push(new PMPoint(-1, 0, 0));
       allpoints.push(new PMPoint(-1, 0, 1));
       allpoints.push(new PMPoint(-1, 1, -1));
       allpoints.push(new PMPoint(-1, 1, 0));
       allpoints.push(new PMPoint(-1, 1, 1));
       allpoints.push(new PMPoint(0, -1, -1));
       allpoints.push(new PMPoint(0, -1, 0));
       allpoints.push(new PMPoint(0, -1, 1));
       allpoints.push(new PMPoint(0, 0, -1));
       allpoints.push(new PMPoint(0, 0, 1));
       allpoints.push(new PMPoint(0, 1, -1));
       allpoints.push(new PMPoint(0, 1, 0));
       allpoints.push(new PMPoint(0, 1, 1));
       allpoints.push(new PMPoint(1, 0, -1));
       allpoints.push(new PMPoint(1, 0, 0));
       allpoints.push(new PMPoint(1, 0, 1));
       allpoints.push(new PMPoint(1, 1, 0));
       allpoints.push(new PMPoint(1, 1, 1));
    
       <!-- Vertex style -->
       var materials = [
          new THREE.MeshBasicMaterial({ color: 0x1EFA1E, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
          new THREE.MeshBasicMaterial({ color: 0x469646, }),
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
       allpoints[7].makesphere(0.02,materials[7]);
       allpoints[8].makesphere(0.02,materials[8]);
       allpoints[9].makesphere(0.02,materials[9]);
       allpoints[10].makesphere(0.02,materials[10]);
       allpoints[11].makesphere(0.02,materials[11]);
       allpoints[12].makesphere(0.02,materials[12]);
       allpoints[13].makesphere(0.02,materials[13]);
       allpoints[14].makesphere(0.02,materials[14]);
       allpoints[15].makesphere(0.02,materials[15]);
       allpoints[16].makesphere(0.02,materials[16]);
       allpoints[17].makesphere(0.02,materials[17]);
       allpoints[18].makesphere(0.02,materials[18]);
       allpoints[19].makesphere(0.02,materials[19]);
       allpoints[20].makesphere(0.02,materials[20]);
    
       for (index = 0; index < allpoints.length; ++index) {
          allpoints[index].add(obj);
       }
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



In order to obtain the integer hull we simply define a new polytope
``$integer_hull`` as the convex hull of all ``LATTICE_POINTS`` contained
in ``$scaled_rs``.

Note that if we give ``POINTS`` (in contrast to ``VERTICES``)
``polymake`` constructs a polytope that is the convex hull of the given
points regardless of whether they are vertices or not. I.e., redundacies
are allowed here.

If you specify ``VERTICES`` you have to make sure yourself that your
points are actually vertices since ``polymake`` does not check this. You
also need to specify the ``LINEALITY_SPACE``, see `Tutorial on
polytopes <apps_polytope>`__.

Linear Programming
~~~~~~~~~~~~~~~~~~

Now that we have constructed a nice integral polytope we want to apply
some linear program to it.

First we define a ``LinearProgram`` with our favourite
``LINEAR_OBJECTIVE``. The linear objective is an given as a vector of
length d+1, d being the dimension of the space. The vector [c0,c1, …,
cd] corresponds to the linear objective c0 + c1x1 + … + cdxd.


::

    polymake> $objective=new LinearProgram<Rational>(LINEAR_OBJECTIVE=>[0,1,1,1]);

Then we define a new polytope, which is a copy of our old one
(``$inter_hull``) with the LP as an additional property.


::

    polymake> $ilp=new Polytope<Rational>(VERTICES=>$integer_hull->VERTICES, LP=>$objective);

|{{ :tutorial:ilp:ilp_min_face.png?200|}}| |{{
:tutorial:ilp:ilp_max_face.png?200|}}|

And now we can perform some computations:

.. |{{ :tutorial:ilp:ilp_min_face.png?200|}}| image:: attachment:ilp_min_face.png
.. |{{ :tutorial:ilp:ilp_max_face.png?200|}}| image:: attachment:ilp_max_face.png


::

    polymake> print $ilp->LP->MAXIMAL_VALUE;
    3




::

    polymake> print $ilp->LP->MAXIMAL_FACE;
    {11}




::

    polymake> $ilp->VISUAL->MIN_MAX_FACE;


.. raw:: html

    <!--
    polymake for andrew
    Thu Mar 28 11:23:53 2019
    ilp
    -->
    
    
    <html>
       <head>
          <title>ilp</title>
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
    <div id="model78121872278"></div>
    
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
       var container = document.getElementById( 'model78121872278' );
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
       var target = document.querySelector('#model78121872278');
    
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
    
       var objectnames = ["ilp"];
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(-1, -1, 0));
       allpoints.push(new PMPoint(-1, 0, -1));
       allpoints.push(new PMPoint(-1, 0, 1));
       allpoints.push(new PMPoint(-1, 1, -1));
       allpoints.push(new PMPoint(-1, 1, 1));
       allpoints.push(new PMPoint(0, -1, -1));
       allpoints.push(new PMPoint(0, -1, 1));
       allpoints.push(new PMPoint(0, 1, -1));
       allpoints.push(new PMPoint(1, 0, -1));
       allpoints.push(new PMPoint(1, 0, 1));
       allpoints.push(new PMPoint(1, 1, 0));
       allpoints.push(new PMPoint(1, 1, 1));
    
       <!-- Vertex style -->
       var materials = [
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
          new THREE.MeshBasicMaterial({ color: 0xFF0000, }),
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
       allpoints[7].makesphere(0.02,materials[7]);
       allpoints[8].makesphere(0.02,materials[8]);
       allpoints[9].makesphere(0.02,materials[9]);
       allpoints[10].makesphere(0.02,materials[10]);
       allpoints[11].makesphere(0.02,materials[11]);
       allpoints[0].makelabel("0");
       allpoints[1].makelabel("1");
       allpoints[2].makelabel("2");
       allpoints[3].makelabel("3");
       allpoints[4].makelabel("4");
       allpoints[5].makelabel("5");
       allpoints[6].makelabel("6");
       allpoints[7].makelabel("7");
       allpoints[8].makelabel("8");
       allpoints[9].makelabel("9");
       allpoints[10].makelabel("10");
       allpoints[11].makelabel("11");
    
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
       faces.vertices.push(allpoints[7].vector);
       faces.vertices.push(allpoints[8].vector);
       faces.vertices.push(allpoints[9].vector);
       faces.vertices.push(allpoints[10].vector);
       faces.vertices.push(allpoints[11].vector);
    
       centroids.push(computeCentroid(faces));
    
       <!-- Facet style -->
       var materials = [
          new THREE.MeshBasicMaterial({ transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, color: 0x77EC9E, }),
          new THREE.MeshBasicMaterial({ transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, color: 0x77EC9E, }),
          new THREE.MeshBasicMaterial({ transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, color: 0x77EC9E, }),
          new THREE.MeshBasicMaterial({ transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, color: 0x77EC9E, }),
          new THREE.MeshBasicMaterial({ transparent: true, opacity: 1, side: THREE.DoubleSide , depthWrite: true, depthTest: true, color: 0xFFFF00, }),
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
       faces.faces.push(new THREE.Face3(8, 7, 10, undefined, undefined, 0));
    
       faces.faces.push(new THREE.Face3(3, 7, 8, undefined, undefined, 1));
       faces.faces.push(new THREE.Face3(3, 8, 5, undefined, undefined, 1));
       faces.faces.push(new THREE.Face3(3, 5, 1, undefined, undefined, 1));
    
       faces.faces.push(new THREE.Face3(5, 6, 0, undefined, undefined, 2));
    
       faces.faces.push(new THREE.Face3(2, 0, 6, undefined, undefined, 3));
    
       faces.faces.push(new THREE.Face3(1, 5, 0, undefined, undefined, 4));
    
       faces.faces.push(new THREE.Face3(4, 3, 1, undefined, undefined, 5));
       faces.faces.push(new THREE.Face3(4, 1, 0, undefined, undefined, 5));
       faces.faces.push(new THREE.Face3(4, 0, 2, undefined, undefined, 5));
    
       faces.faces.push(new THREE.Face3(8, 9, 6, undefined, undefined, 6));
       faces.faces.push(new THREE.Face3(8, 6, 5, undefined, undefined, 6));
    
       faces.faces.push(new THREE.Face3(11, 4, 2, undefined, undefined, 7));
       faces.faces.push(new THREE.Face3(11, 2, 6, undefined, undefined, 7));
       faces.faces.push(new THREE.Face3(11, 6, 9, undefined, undefined, 7));
    
       faces.faces.push(new THREE.Face3(8, 10, 11, undefined, undefined, 8));
       faces.faces.push(new THREE.Face3(8, 11, 9, undefined, undefined, 8));
    
       faces.faces.push(new THREE.Face3(10, 7, 3, undefined, undefined, 9));
       faces.faces.push(new THREE.Face3(10, 3, 4, undefined, undefined, 9));
       faces.faces.push(new THREE.Face3(10, 4, 11, undefined, undefined, 9));
    
    
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
       line.vertices.push(allpoints[8].vector);
       line.vertices.push(allpoints[7].vector);
       line.vertices.push(allpoints[10].vector);
       line.vertices.push(allpoints[8].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[7].vector);
       line.vertices.push(allpoints[8].vector);
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[3].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[5].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[2].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[1].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[1].vector);
       line.vertices.push(allpoints[0].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[4].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[8].vector);
       line.vertices.push(allpoints[9].vector);
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[5].vector);
       line.vertices.push(allpoints[8].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[11].vector);
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[2].vector);
       line.vertices.push(allpoints[6].vector);
       line.vertices.push(allpoints[9].vector);
       line.vertices.push(allpoints[11].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[8].vector);
       line.vertices.push(allpoints[10].vector);
       line.vertices.push(allpoints[11].vector);
       line.vertices.push(allpoints[9].vector);
       line.vertices.push(allpoints[8].vector);
       obj.add(new THREE.Line(line, line_material));
    
       var line = new THREE.Geometry();
       line.vertices.push(allpoints[10].vector);
       line.vertices.push(allpoints[7].vector);
       line.vertices.push(allpoints[3].vector);
       line.vertices.push(allpoints[4].vector);
       line.vertices.push(allpoints[11].vector);
       line.vertices.push(allpoints[10].vector);
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



Hence the LP attains its maximal value 2 on the 2-face spanned by the
vertices 6, 9 and 10.

``polymake`` can visualise the polytope and highlight both its maximal
and minimal face in a different (by default admittedly almost painful
;-) ) colour. Here you see the maximal face ``{6 9 10}`` in red and the
minimal face ``{0 3}`` (on the opposite side of the polytope) in yellow.

Note though that since we started out with a random polytope these
results may vary if we perform the same computations another time on a
different random polytope.


::

    polymake> print $ilp->VERTICES;
    1 -1 -1 0
    1 -1 0 -1
    1 -1 0 1
    1 -1 1 -1
    1 -1 1 1
    1 0 -1 -1
    1 0 -1 1
    1 0 1 -1
    1 1 0 -1
    1 1 0 1
    1 1 1 0
    1 1 1 1





Hilbert bases
~~~~~~~~~~~~~

Finally, we can have ``polymake`` compute and print a Hilbert basis for
the cone spanned by ``$ilp``. Notice that this requires normaliz or 4ti2
to be installed in order to work.


::

    polymake> print $ilp->HILBERT_BASIS;
    1 -1 -1 0
    1 -1 0 -1
    1 -1 0 0
    1 -1 0 1
    1 -1 1 -1
    1 -1 1 0
    1 -1 1 1
    1 0 -1 -1
    1 0 -1 0
    1 0 -1 1
    1 0 0 -1
    1 0 0 0
    1 0 0 1
    1 0 1 -1
    1 0 1 0
    1 0 1 1
    1 1 0 -1
    1 1 0 0
    1 1 0 1
    1 1 1 0
    1 1 1 1





.. raw:: html

    <details><summary><pre style="display:inline"><small>Click here for additional output</small></pre></summary>
    <pre>
    polymake: used package libnormaliz
      Normaliz is a tool for computations in affine monoids, vector configurations, lattice polytopes, and rational cones.
      Copyright by Winfried Bruns, Bogdan Ichim, Christof Soeger.
      http://www.math.uos.de/normaliz/
    
    </pre>
    </details>





