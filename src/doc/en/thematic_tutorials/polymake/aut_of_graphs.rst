.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Automorphisms of Graphs
=======================

First Example: The graph of a square
------------------------------------

Let’s look at the graph of a square. Since a square is a 2-cube, we can
create the polytope and look at its graph:


.. link

.. CODE-BLOCK:: perl

    polymake> $c=cube(2);
    polymake> $c->GRAPH->VISUAL;


.. raw:: html

    <!--
    polymake for andrew
    Thu Apr  4 11:21:35 2019
    unnamed
    -->
    
    
    <html>
       <head>
          <title>unnamed</title>
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
    <div id="model7073906164"></div>
    
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
       var container = document.getElementById( 'model7073906164' );
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
       var target = document.querySelector('#model7073906164');
    
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
    
       var objectnames = ["unnamed__1"];
       var obj = new THREE.Object3D();
       var allpoints = [];
       allpoints.push(new PMPoint(-0.326103, -1.28762, -0.341913));
       allpoints.push(new PMPoint(0.365144, -0.424679, 1.25201));
       allpoints.push(new PMPoint(-0.365144, 0.424679, -1.25201));
       allpoints.push(new PMPoint(0.326103, 1.28762, 0.341913));
    
       <!-- Vertex style -->
       var points_material = new THREE.MeshBasicMaterial ( {color: 0xFF0000, } );
    
       points_material.side = THREE.DoubleSide;
       points_material.transparent = true;
    
       <!-- POINTS -->
       allpoints[0].makesphere(0.02,points_material);
       allpoints[1].makesphere(0.02,points_material);
       allpoints[2].makesphere(0.02,points_material);
       allpoints[3].makesphere(0.02,points_material);
       allpoints[0].makelabel("0");
       allpoints[1].makelabel("1");
       allpoints[2].makelabel("2");
       allpoints[3].makelabel("3");
    
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




.. figure:: attachment:square.png
   :alt: {{:tutorial:square.png?200|}}

   {{:tutorial:square.png?200|}}

To study the automorphisms of this graph, we create a ``props::Graph``
object refering to the ``C++`` class named ``Graph`` (see the `tutorial
on graphs <apps_graph>`__ for more details):


.. link

.. CODE-BLOCK:: perl

    polymake> $g=new props::Graph($c->GRAPH->ADJACENCY);

The picture of the graph shows that the node with label 0 is adjacent to
the nodes 1 and 2, Node 1 is adjacent to 0 and 3, and so on. For the
complete adjacency information you can print ``$c->GRAPH->ADJACENCY`` or
just the props::Graph object ``$g``:


.. link

.. CODE-BLOCK:: perl

    polymake> print rows_labeled($g);
    0:1 2
    1:0 3
    2:0 3
    3:1 2





Now, we compute the generators of the automorphism group of this graph
(see the `tutorial on groups <apps_group>`__ for more info):


.. link

.. CODE-BLOCK:: perl

    polymake> $aut=automorphisms($g);


.. raw:: html

    <details><summary><pre style="display:inline"><small>Click here for additional output</small></pre></summary>
    <pre>
    polymake: used package nauty
      Computation of automorphism groups of graphs.
      Copyright by Brendan McKay and Adolfo Piperno.
      http://pallini.di.uniroma1.it/
    </pre>
    </details>




In this case, the automorphism group has two generators:


.. link

.. CODE-BLOCK:: perl

    polymake> print $aut;
    0 2 1 3
    1 0 3 2





Each generator is a permutation on the nodes. The first generator fixes
the nodes 0 and 3, and exchanges the nodes 1 and 2, i.e., it describes
the reflection along the diagonal through 0 and 3. The second generator
is the reflection along the horizontal line.

In order to be able to work with the group, we create a new Group
object, which lives in the application ``group``:


.. link

.. CODE-BLOCK:: perl

    polymake> $action = new group::PermutationAction(GENERATORS => $aut);
    polymake> $autgroup = new group::Group(PERMUTATION_ACTION => $action);

Now we can ask for basic properties of the group, e.g., the number of
elements:


.. link

.. CODE-BLOCK:: perl

    polymake> print $autgroup->ORDER;
    8




.. raw:: html

    <details><summary><pre style="display:inline"><small>Click here for additional output</small></pre></summary>
    <pre>
    polymake: used package permlib
      A callable C++ library for permutation computations. 
      Written by Thomas Rehn.
      http://www.math.uni-rostock.de/~rehn/software/permlib.html 
    
    </pre>
    </details>




Sometimes, it is useful to know which elements of the group fix a
specific set of indices, that is, we are interested in the subgroup
which is the stabilizer of the given set. In the first case, we just fix
the index 0:


.. link

.. CODE-BLOCK:: perl

    polymake> $s0=new Set<Int>(0);
    polymake> $stab0=group::stabilizer_of_set($action,$s0);

We learn that the node 0 is only fixed by the permutation ``0 2 1 3``:


.. link

.. CODE-BLOCK:: perl

    polymake> print $stab0->ORDER;
    2




.. link

.. CODE-BLOCK:: perl

    polymake> print $stab0->PERMUTATION_ACTION->GENERATORS;
    0 2 1 3





In the second case, we look at the subgroup which leaves the set
``{1,2}`` invariant:


.. link

.. CODE-BLOCK:: perl

    polymake> $s12=new Set<Int>(1,2);
    polymake> $stab12=group::stabilizer_of_set($action,$s12);

Now, we obtain a group of order 4:


.. link

.. CODE-BLOCK:: perl

    polymake> print $stab12->ORDER;
    4




.. link

.. CODE-BLOCK:: perl

    polymake> print $stab12->PERMUTATION_ACTION->GENERATORS;
    3 1 2 0
    0 2 1 3





Finally, we compute the orbits of the indices under the three different
groups:


.. link

.. CODE-BLOCK:: perl

    polymake> print $stab0->PERMUTATION_ACTION->ORBITS;
    {0}
    {1 2}
    {3}





.. link

.. CODE-BLOCK:: perl

    polymake> print $stab12->PERMUTATION_ACTION->ORBITS;
    {0 3}
    {1 2}





.. link

.. CODE-BLOCK:: perl

    polymake> print $autgroup->PERMUTATION_ACTION->ORBITS;
    {0 1 2 3}



