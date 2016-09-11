var Field = function(N) {
  this.N = N;
  this.values = new Array(N+2);
  for(var i=0; i<this.values.length; i++) {
    this.values[i] = new Array(N+2);
    for (var j=0; j<N+2; j++) {
      this.values[i][j] = 0.0;
    }
  }
};

Field.prototype.swap = function(other) {
  var tmp = other.values;
  other.values = this.values;
  this.values = tmp;
};

Field.prototype.clear = function() {
  for (var i=0; i<this.values.length; i++) {
    for (var j=0; j<this.values.length; j++) {
      this.values[i][j] = 0;
    }
  }

};

Field.prototype.addSource = function(source, dt) {
  for (var i=0; i<this.values.length; i++) {
    for (var j=0; j<this.values.length; j++) {
      this.values[i][j] += dt*source.values[i][j];
    }
  }
};

Field.prototype.setBoundary = function(b) {
  var N = this.N;

  for (var i=1; i<=N; i++) {
    this.values[0][i]   = (b == 1 ? -1.0: 1.0) * this.values[1][i];
    this.values[N+1][i] = (b == 1 ? -1.0: 1.0) * this.values[N][i];
    this.values[i][0]   = (b == 2 ? -1.0: 1.0) * this.values[i][1];
    this.values[i][N+1] = (b == 2 ? -1.0: 1.0) * this.values[i][N];
  }

  this.values[0][0]     = 0.5 * (this.values[1][ 0 ] + this.values[ 0 ][1]);
  this.values[0][N+1]   = 0.5 * (this.values[1][N+1] + this.values[ 0 ][N]);
  this.values[N+1][0]   = 0.5 * (this.values[N][ 0 ] + this.values[N+1][1]);
  this.values[N+1][N+1] = 0.5 * (this.values[N][N+1] + this.values[N+1][N]);
};

var Solver = function(N, diff, visc) {
  this.N = N;
  this.diff = diff || 0.0;
  this.visc = visc || 0.0;
  this.u  = new Field(this.N);
  this.u_prev = new Field(this.N);
  this.v  = new Field(this.N);
  this.v_prev = new Field(this.N);
  this.density = new Field(this.N);
  this.density_prev = new Field(this.N);
};

Solver.prototype.reset = function() {
  this.u.clear();
  this.u_prev.clear();
  this.v.clear();
  this.v_prev.clear();
  this.density.clear();
  this.density_prev.clear();
};

Solver.prototype.inject = function(x, y, density, force, dx, dy) {
  this.u_prev.values[x][y] = force * dx;
  this.v_prev.values[x][y] = force * dy;
  this.density_prev.values[x][y] = density;
}

Solver.prototype.linearSolve = function(x, x0, b, a, c) {
  for (var k=0; k<20; k++) {
    for (var i=1; i<=this.N; i++) {
      for (var j=1; j<=this.N; j++) {
        x.values[i][j] = (x0.values[i][j] +
                          a*(x.values[i-1][j] + x.values[i+1][j] +
                             x.values[i][j-1] + x.values[i][j+1]))/c;
      }
    }
    x.setBoundary(b);
  }
};

Solver.prototype.diffuse = function(x, x0, b, diff, dt) {
  var a = dt*diff*this.N*this.N;
  this.linearSolve(x, x0, b, a, 1+4*a);
};

Solver.prototype.advect = function(d, d0, u, v, b, dt) {
  var i0, j0, i1, j1;
  var x, y, s0, t0, s1, t1, dt0;
  var N = this.N;

  dt0 = dt*N;

  for(var i=1; i<=this.N; i++) {
    for(var j=1; j<=this.N; j++) {
      x = i - dt0*u.values[i][j];
      y = j - dt0*v.values[i][j];
      if (x<0.5) x=0.5;
      if (x>N+0.5) x=N+0.5;
      i0 = Math.floor(x); i1 = i0+1;

      if (y<0.5) y=0.5;
      if (y>(N+0.5)) y = N+0.5;
      j0 = Math.floor(y); j1 = j0+1;

      s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
      d.values[i][j] = s0*(t0*d0.values[i0][j0] + t1*d0.values[i0][j1]) +
        s1*(t0*d0.values[i1][j0] + t1*d0.values[i1][j1]);
    }
  }

  d.setBoundary(b);
};

Solver.prototype.project = function(u, v, p, div) {
  var i, j;

  for(i=1; i<=this.N; i++) {
    for(j=1; j<=this.N; j++) {
      div.values[i][j] = -0.5*(u.values[i+1][j]-u.values[i-1][j]+v.values[i][j+1]-v.values[i][j-1])/this.N;
      p.values[i][j] = 0;
    }
  }

  div.setBoundary(0);
  p.setBoundary(0);

  this.linearSolve(p, div, 0, 1, 4);

  for(i=1; i<=this.N; i++) {
    for(j=1; j<=this.N; j++) {
      u.values[i][j] -= 0.5*this.N*(p.values[i+1][j]-p.values[i-1][j]);
      v.values[i][j] -= 0.5*this.N*(p.values[i][j+1]-p.values[i][j-1]);
    }
  }

  u.setBoundary(1);
  v.setBoundary(2);
};

Solver.prototype.densityStep = function(x, x0, u, v, diff, dt) {
  x.addSource(x0, dt);
  x0.swap(x);
  this.diffuse(x, x0, 0, diff, dt);
  x0.swap(x);
  this.advect(x, x0, u, v, 0, dt);
};

Solver.prototype.velocityStep = function(u, v, u0, v0, visc, dt) {
  u.addSource(u0, dt);
  v.addSource(v0, dt);
  u0.swap(u);
  this.diffuse(u, u0, 3, visc, dt);
  v0.swap(v);
  this.diffuse(v, v0, 3, visc, dt);
  this.project(u, v, u0, v0);
  u0.swap(u);
  v0.swap(v);
  this.advect(u, u0, u0, v0, 1, dt);
  this.advect(v, v0, u0, v0, 2, dt);
  this.project(u, v, u0, v0);
};

Solver.prototype.step = function(dt) {
  dt = dt || 0.05;
  this.velocityStep(this.u, this.v, this.u_prev, this.v_prev, this.visc, dt);
  this.densityStep(this.density, this.density_prev, this.u, this.v, this.diff, dt);

  // density_prev must be cleared between runs
  this.density_prev.clear();
  this.u_prev.clear();
  this.v_prev.clear();
};

var solver;
var container;
var camera, cameraRTT, scene, sceneRTT, renderer, context;
var renderTexture;
var N = 128; // resolution of sim
var mouseX = null, mouseY = null;
var dataMaterial, quad, dataTexture;
var dt = 0.01;
var imageData = new Float32Array(N*N*3);
var vizColorMap;
var densityColorMap;

function lerp(a, b, t) {
  return (1-t)*a+t*b;
}

var ColorMap = function() {
  this.stops = [];
};

ColorMap.prototype.addColorStop = function(pos, r, g, b) {
  for(var i=0; i<this.stops.length; i++) {
    if (this.stops[i].pos > pos) {
      i++;
      break;
    }
  }

  this.stops[i] = {
    pos: pos,
    r: r,
    g: g,
    b: b
  };
}

ColorMap.prototype.sample = function(t) {
  var a=0, b=0;
  t = Math.min(1, Math.max(0, t));

  for(var i=0; i<this.stops.length; i++) {
    if (t < this.stops[i].pos) {
      a=i>0? i-1 : 0;
      b=i;
      break;
    }
  }

  if (b == 0) {
    a = b = this.stops.length - 1;
  }

  if (a==b) {
    return [this.stops[a].r, this.stops[a].g, this.stops[a].b];
  }

  var A=this.stops[a];
  var B=this.stops[b];
  var t1 = (t-A.pos)/(B.pos-A.pos);
  return [lerp(A.r, B.r, t1), lerp(A.g, B.g, t1), lerp(A.b, B.b, t1)];
};

(function main() {
  THREE.ImageUtils.crossOrigin = "anonymous";
  init();
  initData(solver);
  animate();
})();


function handleKeyboardInput(e) {
  switch(e.keyCode) {
    case 114:
      initData(solver);
      break;
    default:
      break;
  }
}

function initData(solver) {
  solver.reset();

  var size = (solver.N+2)^2;
  for(var i=0; i<solver.N+2; i++) {
    for(var j=0; j<solver.N+2; j++) {
      var t = 1.0 - j/solver.N;
      solver.density.values[i][j] = densityColorMap.sample(t)[0];
    }
  }
}

function init() {
  solver = new Solver(N, 0.0, 0.0005);

  container = document.getElementById("container");
  camera = new THREE.PerspectiveCamera(30, window.innerWidth / window.innerHeight, 1, 10);
  camera.position.z = 2;

  scene = new THREE.Scene();

  cameraRTT = new THREE.OrthographicCamera(-0.5, 0.5, 0.5, -0.5, 1, 10);
  cameraRTT.position.z = 2;
  sceneRTT = new THREE.Scene();

  renderTexture = new THREE.WebGLRenderTarget(512, 512, { minFilter: THREE.LinearFilter, magFilter: THREE.NearestFilter, format: THREE.RGBFormat });

  dataTexture = new THREE.DataTexture(imageData, N, N, THREE.RGBFormat, THREE.FloatType);
  dataTexture.needsUpdate = true;

  dataMaterial = new THREE.MeshBasicMaterial({ map: dataTexture });

  var plane = new THREE.PlaneBufferGeometry(1, 1, 1, 1); // TODO: make this resize on screen resize

  quad = new THREE.Mesh(plane, dataMaterial);
  sceneRTT.add(quad);

  // Background
  {
    var bgTex = THREE.ImageUtils.loadTexture("https://drive.google.com/file/d/0BwpF5xsAFjCYSzV5dHdWbVc1UUk/view?usp=sharing");
    var aspect = 1.6/2.5;
    var scale = 2.5;
    var bgMaterial = new THREE.MeshBasicMaterial({ map: bgTex, depthWrite: false });
    var bgGeo = new THREE.PlaneBufferGeometry(scale, scale*aspect, 1, 1);
    var bgMesh = new THREE.Mesh(bgGeo, bgMaterial);
    scene.add(bgMesh);
  }

  // Glow
  // {
  //   var size = 2.0;
  //   var glowMat = new THREE.ShaderMaterial({
  //     vertexShader: document.getElementById("shared_vert").textContent,
  //     fragmentShader: document.getElementById("glow_frag").textContent,
  //     blending: THREE.AdditiveBlending,
  //     depthWrite: false,
  //     transparent: true
  //   });
  //   var glowGeo = new THREE.PlaneBufferGeometry(size, size, 1, 1);
  //   var glowMesh = new THREE.Mesh(glowGeo, glowMat);
  //   scene.add(glowMesh);
  // }

  // Tea
  {
    var teaMat = new THREE.ShaderMaterial({
      vertexShader: document.getElementById("shared_vert").textContent,
      fragmentShader: document.getElementById("fluid_frag").textContent,
      uniforms: {
        "texture": {
          type: "t",
          value: renderTexture
        },
        "reflection": {
          type: "t",
          value: THREE.ImageUtils.loadTexture("https://dl.dropboxusercontent.com/u/8604128/tea/rainbow-nebula-sphere-blur.jpg")
        }
      }
    });

    var teaGeo = new THREE.CylinderGeometry(1, 1, 1, 30, 10, true);
    var teaMesh = new THREE.Mesh(teaGeo, teaMat);
    teaMesh.rotation.y = 180;
    scene.add(teaMesh);
  }

  // Cup
 // {
 //    var cupMaterial = new THREE.ShaderMaterial({
 //      vertexShader: document.getElementById("shared_vert").textContent,
 //      fragmentShader: document.getElementById("cup_frag").textContent,
 //      blending: THREE.AdditiveBlending,
 //      transparent: true,
 //      depthWrite: false,
 //      uniforms: {
 //        "texture": {
 //          type: "t",
 //          value: THREE.ImageUtils.loadTexture("https://dl.dropboxusercontent.com/u/8604128/tea/rainbow-nebula-sphere.jpg")
 //        }
 //      }
 //    });
 //    var cupGeo = new THREE.CylinderGeometry(0.32, 0.22, 0.82, 32, 2, true);
 //    var cupMesh = new THREE.Mesh(cupGeo, cupMaterial);
 //    scene.add(cupMesh);
 //  }

  // Scene renderer
  renderer = new THREE.WebGLRenderer();
  renderer.setSize(window.innerWidth, window.innerHeight);
  renderer.autoClear = false;
  container.appendChild(renderer.domElement);

  // Create color map for the tea/milk gradient
  vizColorMap = new ColorMap();
  vizColorMap.addColorStop(0.3, 0, 0, 0);
  vizColorMap.addColorStop(0.3, 0 ,0, 238);
  vizColorMap.addColorStop(0.3, 0 ,0, 255);
  vizColorMap.addColorStop(1.0, 1.0, 1.0, 1.0);

  // Create color map defining initial data breakdown
  densityColorMap = new ColorMap();
  densityColorMap.addColorStop(0.0, 0, 0, 0);
  densityColorMap.addColorStop(0.8, 0.5, 0.5, 0.5);
  densityColorMap.addColorStop(0.803, 1.0, 1.0, 1.0);

  // Events
  document.addEventListener('mousemove', onMouseMove, false);
  document.addEventListener('keypress', handleKeyboardInput);
  window.addEventListener('resize', resize, false);
}

Number.prototype.clamp = function(low, high) {
  return Math.min(Math.max(low, this), high);
}

Number.prototype.lerp = function(low, high) {
  return (1.0-this)*low + this * high;
}

function onMouseMove(event) {
  var bbox = container.getBoundingClientRect();
  var newMouseX = (event.clientX - bbox.left)/(bbox.right - bbox.left);
  var newMouseY = (event.clientY - bbox.top)/(bbox.bottom - bbox.top);

  // remap to fit cup area
  var startX = 0.33;
  var endX = 0.66;
  newMouseX = ((newMouseX.clamp(startX, endX) - startX)/(endX-startX)).lerp(0.22, 0.52);


  mouseX = mouseX || newMouseX;
  mouseY = mouseY || newMouseY;

  var dx = newMouseX - mouseX;
  var dy = newMouseY - mouseY;

  mouseX = newMouseX;
  mouseY = newMouseY;

  var x = Math.floor(mouseX * solver.N);
  var y = Math.floor(mouseY * solver.N);
  solver.inject(x, y, 0.0, 50.0, dx*solver.N, dy*solver.N);
}

var lastTime = Date.now();
var startTime = Date.now();
var elapsedTime = 0;
function animate() {
  requestAnimationFrame(animate);

  var now = Date.now();
  var dt = (now - lastTime)/1000.0;

  lastTime = now;
  elapsedTime += dt;

  solver.step(dt);

  updateTexture(solver, dataTexture, vizColorMap);
  render();
}

function updateTexture(solver, texture, colorMap) {
  for(var i=1; i<=solver.N; i++) {
    for(var j=1; j<=solver.N; j++) {
      var x = i-1;
      var y = j-1;
      var color = colorMap.sample(solver.density.values[i][j]);
      imageData[(x + solver.N * y)*3 + 0] = color[0];
      imageData[(x + solver.N * y)*3 + 1] = color[1];
      imageData[(x + solver.N * y)*3 + 2] = color[2];
    }
  }

  texture.needsUpdate = true;
}

function render() {
  var time = Date.now();
  renderer.clear();

  // render fluid simulation to texture:
  renderer.render(sceneRTT, cameraRTT, renderTexture, true);

  // render scene:
  renderer.render(scene, camera);
}

function resize() {
  var width = window.innerWidth;
  var height = window.innerHeight;
  var aspect = width/height;
  camera.aspect = aspect;
  camera.updateProjectionMatrix();
  renderer.setSize(width, height);
}
