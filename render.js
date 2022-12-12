// values declared
let scene, camera, renderer, cube;
// documentation for reference in index.html
let kikuchi = document.getElementById("ThreeDiv");

// defining constants

// default quaternion values
var quaternion = new THREE.Quaternion();
var q1 = 0;
var q2 = 0;
var q3 = 0;
var q4 = 1;

//unit vectors
const x = new THREE.Vector3( 1, 0, 0 );
const y = new THREE.Vector3( 0, 1, 0 );
const z = new THREE.Vector3( 0, 0, 1 );

//const radius = 0.3; //radius of the cylinders

// lattice vectors
const a = 3.19e-10;
const c = 5.18e-10;
const u = (3 * c * c + 4 * a * a)/(12 * c * c);

//const W = (3*(a**2))/(4*(c**2));

const n = 1; // order of diffraction
const lightspeed = 299792458; //speed of light
const me = 9.1094897*10**-31; //mass of an electron
const e = 1.60217733*10**-19; //charge of an electron
const plancksconstant = 6.62607004*10**-34; //plancks constant
const a0 = 5.29177210903*10**-11; //Bohr radius
var beamV = 20000; // accelerating voltage of electron beam

//coordinates of Ga atoms in units cell:
var Ga1x = 1 / 3;
var Ga1y = 2 / 3;
var Ga1z = 0;
var Ga2x = 2 / 3;
var Ga2y = 1 / 3;
var Ga2z = 1 / 2;
//coordinates of N atoms (offset from Ga atoms by u in c direction):
var N1x = 1 / 3;
var N1y = 2 / 3;
var N1z = u;
var N2x = 2 / 3;
var N2y = 1 / 3;
var N2z = 1 / 2 + u;

const reflectors = [];

function loopHKL(maxHKL) {
  
  //Loop over h k l values such that none of h k i(=0-h-k) or l has an absolute value greater than maxHKL.
  //Keep h >= k >= i, so no need to start at h = maxHKL (which will be covered by the i index).
  //Instead start at h = max/HKL (rounded up):
  maxh = Math.round(maxHKL / 2)
  for (let h = maxh; h >= 0; h--){
    for (let k = h; k >= 0; k--){
      let i = -h - k;
      //Avoid i going to -(maxHKL+1) when maxHKL is odd:
      if (-i <= maxHKL) {
        for (let l = maxHKL; l >= -maxHKL; l--){
          //Don't include (0 0 0):
          if (h!=0 || k!=0 || l!=0){ 
            //console.log("(" + h + " " + k + " " + i + " " + l + ")" );
            let twod = calc2d(h,k,l);
            //console.log("2d=" + twod);
            let fGa = eFormFac(31, twod);
            //console.log("fGa=" + fGa);
            let fN = eFormFac(7, twod);
            //console.log("fN=" + fN);
            let SF = calcSF(h, k, l, twod);
            //console.log("StrucFac=" + SF);
            let mult = multiplicity(h, k, l);
            //console.log("Mult=" + mult);
            let reflector = {h:h, k:k, l:l, F:SF, mult:mult};
            reflectors.push(reflector); 
          }
        }
      }
    }
  }
  reflectors.sort(function(a, b){return b.F - a.F});
  reflectors.forEach(r => {
    let i = 0 - r.h - r.k;
    let rtext = ("(" + r.h + " " + r.k + " " + i + " " + r.l + ")\t|F|=" + r.F + "\tMult=" + r.mult + "\n");
    //console.log(rtext);
  });
  console.log(reflectors);
}

function calcSF(h, k, l, twod){
  let X = 4 /3 * h + 2 / 3 * k + l;
  //Check if this reflection is forbidden:
  if (Math.round((X + 1) / 2) == (X + 1) / 2){
    //console.log("Forbidden");
    return 0;
  } else {
    //console.log("Allowed");
    let eFFGa = eFormFac(31, twod);
    let eFFN = eFormFac(7, twod);
    let real = eFFGa * (1 + Math.cos(Math.PI * X))
      + eFFN * (Math.cos(Math.PI * 2 * u * l) + Math.cos(Math.PI * (X + 2 * u * l)));
    let imag = eFFGa * (Math.sin(-Math.PI * X))
      + eFFN * (Math.sin(-Math.PI * 2 * u * l) + Math.sin(-Math.PI * (X + 2 * u * l)));
    //console.log(real + "+" + imag + "i");
    return Math.sqrt(real * real + imag * imag);
  }
}

function xFormFac(Z, twod){
  //Calculates form factor for atom of atomic number Z
  //q = (sin theta)/lambda = 1/2d
  twod2 = twod *twod;
  if(Z == 31){
    //Ga:
    return 15.2354 * Math.exp(-3.0669 / twod2)
      + 6.7006 * Math.exp(-0.2412 / twod2)
      + 4.3591 * Math.exp(-10.7805 / twod2)
      + 2.9623 * Math.exp(-61.4135 / twod2)
      + 1.7189;
  } else if (Z == 7){
    //N:
    return 12.2126 * Math.exp(-0.00570 / twod2)
      + 3.1322 * Math.exp(-9.8933 / twod2)
      + 2.0125 * Math.exp(-28.9975 / twod2)
      + 1.1663 * Math.exp(-0.5826 / twod2)
      - 11.529;
  }
}

function multiplicity(h, k, l){
  let mult = 12;
  //absolute values of h, k, i:
  let ha = Math.abs(h);
  let ka = Math.abs(k);
  let ia = Math.abs(-h - k);
  
  if (h == 0 && k == 0) return 1;
  if (ha == ka || ka == ia || ia == ha) mult /= 2;
  //if (l == 0) mult /= 2;
  return mult;
  
}

function eFormFac(Z, twod){
  let xFF = xFormFac(Z, twod);
  return 0.02393 * (Z - xFF) * twod *twod;
}

function calc2d(h, k, l){
  return 2 * a / Math.sqrt(4 / 3 * (h * h + h * k + k * k) + Math.pow(a * l / c, 2));
}

// calculating the wavelength of the beam
var lamda = plancksconstant/(Math.sqrt(2*me*e*beamV*(1+(e/2*me*(lightspeed)**2)*beamV)));

function kinematic(){
  loopHKL(6)
  let maximumS = reflectors[0].F
  let rads = 0.301;
  for(j = 0; j<reflectors.length;j++){
    if(reflectors[j].F>0.01*maximumS){
      let intensity = reflectors[j].F/maximumS
      console.log(intensity)
      let h = reflectors[j].h
      let k = reflectors[j].k
      let i = -(h+k)
      let l = reflectors[j].l
      let radius = rads + 0.001*j
      Highlight(h,k,i,l,intensity,radius,reflectors[j].mult)
    }

  }
}

function Highlight(h,k,i,l,intensity,radius,bandnumber) {
  //Creates a band from an input by calculating positon
  let opacitie = Number(intensity)
  let width = bandWidth(h,k,l,radius);
  //let bandnumber = Number(Multiplicity(h,k,i));
  var cylGeometry = new THREE.CylinderGeometry(radius, radius, width, 30, 30, true);
  var cylMaterial = new THREE.MeshBasicMaterial({color: 0xFFFFFF, side: THREE.DoubleSide, transparent: true, opacity: opacitie});
  const N2A = N2Aangle(h,k,i,l)
  const N2C = N2Cangle(h,k,i,l)
  for (let count = 0; count < (bandnumber); count++) {
    cylinder = new THREE.Mesh(cylGeometry, cylMaterial);
    cylinder.rotateOnWorldAxis(x,Math.PI/2);
    cylinder.rotateOnWorldAxis(x,N2C);
    if (count>5){
      cylinder.rotateOnWorldAxis(z,-N2A-Math.PI*2*count/6)
    } else{
      cylinder.rotateOnWorldAxis(z,N2A+Math.PI*2*count/6)
    }
    scene.add(cylinder);
  }
}


function init() {
    // sets scene and camera angle
    scene = new THREE.Scene();
    camera = new THREE.PerspectiveCamera( 90, 400 / 400, 0.0001, 10 );
    const axesHelper = new THREE.AxesHelper( 5 );
    scene.add(axesHelper);
  
    // begins rendering process, edges smoothed
    renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setSize(400, 400);
    kikuchi.appendChild(renderer.domElement);
  
     //cube defined
    const geometry = new THREE.BoxGeometry();
    const loader = new THREE.TextureLoader();
    const materials = [
      new THREE.MeshBasicMaterial({map: loader.load(`./Model/Black.png`), side: THREE.DoubleSide}),
      new THREE.MeshBasicMaterial({map: loader.load(`./Model/Black.png`), side: THREE.DoubleSide}),
      new THREE.MeshBasicMaterial({map: loader.load(`./Model/Black.png`), side: THREE.DoubleSide}),
      new THREE.MeshBasicMaterial({map: loader.load(`./Model/Black.png`), side: THREE.DoubleSide}),
      new THREE.MeshBasicMaterial({map: loader.load(`./Model/Black.png`), side: THREE.DoubleSide}),
      new THREE.MeshBasicMaterial({map: loader.load(`./Model/Black.png`), side: THREE.DoubleSide}),
    ];
    cube = new THREE.Mesh(geometry, materials);
  
    // generates cube
    //scene.add(cube);

    // z-axis position of camera, centre of cube at position 0
    camera.position.z = 0;

}

function TextureChange() {
  scene.remove(cube)
  let CIF = document.getElementById("Cif").value;
  let AcceleratingVoltage = document.getElementById("EAV").value;
  let Texture = `${CIF}/${AcceleratingVoltage}keV`
  const geometry = new THREE.BoxGeometry();
  const loader = new THREE.TextureLoader();
  const materials = [
    new THREE.MeshBasicMaterial({map: loader.load(`./Model/${Texture}/18090n90.png`), side: THREE.DoubleSide}),
    new THREE.MeshBasicMaterial({map: loader.load(`./Model/${Texture}/909090.png`), side: THREE.DoubleSide}),
    new THREE.MeshBasicMaterial({map: loader.load(`./Model/${Texture}/090180.png`), side: THREE.DoubleSide}),
    new THREE.MeshBasicMaterial({map: loader.load(`./Model/${Texture}/n90900.png`), side: THREE.DoubleSide}),
    new THREE.MeshBasicMaterial({map: loader.load(`./Model/${Texture}/18000.png`), side: THREE.DoubleSide}),
    new THREE.MeshBasicMaterial({map: loader.load(`./Model/${Texture}/01800.png`), side: THREE.DoubleSide}),
  ];
  cube = new THREE.Mesh(geometry, materials);
  
  // generates cube
  scene.add(cube);
  Jmol.script(myJmol,`load ./Model/${CIF}.cif {445,665,-1}`)

}

//calculates the width of the kikuchi band through Bragg's Law
function bandWidth(h,k,l,radius) {
  let AcceleratingVoltage = document.getElementById("EAV").value;
  var beamV = AcceleratingVoltage*1000;
  var lamda = plancksconstant/(Math.sqrt(2*me*e*beamV*(1+(e/2*me*(lightspeed)**2)*beamV)));
  let d = a/(Math.sqrt(4/3*((h**2) + h*k + (k**2))+((a**2)/(c**2))*(l**2)))
  let BraggAngle = Math.asin((n*lamda)/(2*(d)));
  let width = 2*radius*Math.tan(BraggAngle);
  return width;
}

function N2Aangle(h,k,i,l) {
  // calculates the angle between the normal of the plane and the direction of the a-axis [2,-1,-1,0]
  W = l*((3*a**2)/(2*c**2))*0 //when this term is multiplied by zero, the program works.
  const num = (a**2)*(3*(2*h-k)+(3/2)*(2*k-h)) // numerator of the equation
  const den = 3*a*Math.sqrt((3*(a**2)*((h**2)+h*k+(k**2)))+((c**2)*(W**2))) // denomenator of equation
  const angle = Math.acos(num/den); //Angle between normal and direction of A-axis [2,-1,-1,0]
  console.log("Normal to A axis:",angle);
  if (isNaN(angle)){
    return 0;
  }
  return angle;
}

function N2Cangle(h,k,i,l) {
  // calculates the angle between the normal of the plane and the direction of the c-axis [0,0,0,1]
  W = ((3*a**2)/(2*c**2))*l
  let angle = Math.acos((W*(c**2))/(c*Math.sqrt(3*(a**2)*((h**2)+h*k+(k**2))+(c**2)*(W**2))));
  console.log("Normal to C axis",angle);
  if (isNaN(angle)){
    return 0;
  }
  return angle;
}

//hammond the basics of crystallography and diffraction appendix 4 - angle between planes equation - suspect problems - didnt work

function HighlightBands(h,k,i,l) {
  //Creates a band from an input by calculating positon
  let width = bandWidth(h,k,l,0.3);
  let bandnumber = Number(Multiplicity(h,k,i));
  var cylGeometry = new THREE.CylinderGeometry(0.3, 0.3, width, 30, 30, true);
  var cylMaterial = new THREE.MeshBasicMaterial({color: 0xFF0000, side: THREE.DoubleSide, transparent: true, opacity: 0.3});
  const N2A = N2Aangle(h,k,i,l)
  const N2C = N2Cangle(h,k,i,l)
  for (let count = 0; count < (bandnumber); count++) {
    cylinder = new THREE.Mesh(cylGeometry, cylMaterial);
    cylinder.rotateOnWorldAxis(x,Math.PI/2);
    cylinder.rotateOnWorldAxis(x,N2C);
    if (count>5){
      cylinder.rotateOnWorldAxis(z,-N2A-Math.PI*2*count/6)
    } else{
      cylinder.rotateOnWorldAxis(z,N2A+Math.PI*2*count/6)
    }
    cylinder.name = count
    scene.add(cylinder);
  }
}

function Multiplicity(h,k,i) {
  // caluclates the multiplicity - the number of bands that relate to a specified plane
  equivplanes = permutator([h,k,i]);
  m = equivplanes.length;
  //inverse of permutations added to same array
  for (let count = 0; count < m; count++){
    let negplane = [];
    plane = equivplanes[count];
    for (let ncount = 0; ncount < 3; ncount++){
      negplane[ncount] = plane[ncount]*-1
    }
    equivplanes.push(negplane)
  }
  // remove duplicates
  equivplanes = eliminateDuplicates(equivplanes)
  console.log(equivplanes)
  m = equivplanes.length
  return m;
}

function eliminateDuplicates(arr) {
  let i,
      len = arr.length,
      out = [],
      obj = {};
  for (i = 0; i < len; i++) {
    obj[arr[i]] = 0;
  }
  for (i in obj) {
    out.push(i);
  }
  return out;
}

function permutator(inputArr) {
  var results = [];
  function permute(arr, memo) {
    var cur, memo = memo || [];
    for (var i = 0; i < arr.length; i++) {
      cur = arr.splice(i, 1);
      if (arr.length === 0) {
        results.push(memo.concat(cur));
      }
      permute(arr.slice(), memo.concat(cur));
      arr.splice(i, 0, cur[0]);
    }
    return results;
  }

  return permute(inputArr);
}

function removebands() {
  //removes cylinders that highlight the kikuchi bands within the 101 plane
  for (let j = 0; j < 30 ;j++) {
    scene.remove(scene.getObjectByName(j));
    scene.remove(scene.getObjectByName(j));
  }
}

function animate() {
    // controls animation of cube
    requestAnimationFrame(animate);
    quaternion.set(q1,q2,q3,q4);
    quaternion.normalize();
    // camera rotation set by quaternion values
    camera.rotation.setFromQuaternion(quaternion);
    renderer.render(scene, camera);
  }

function locate(a,b){
  // relays the rotation information from the unit cell to the pattern by changing the quarternion values
  var str = b;
  if (str.slice(0,1)=="{") {
    var qarray = str.slice(1,str.length-1).split(" ");
    q1 = qarray[0];
    q2 = qarray[1];
    q3 = qarray[2];
    q4 = qarray[3];
  }
  animate();
}

function reset(){
  q1 = 0;
  q2 = 0;
  q3 = 0;
  q4 = 1;
  
  animate();
}


function SubmitVals() {
  Jmol.script(myJmol,'isosurface p4 delete')
  removebands();
  const HInputs = document.getElementById("HInput").value;
  const KInputs = document.getElementById("KInput").value;
  const IInputs = -(Number(HInputs)+Number(KInputs))
  document.getElementById('IInputs').innerHTML = IInputs;
  const LInputs = document.getElementById("LInput").value;
  Jmol.script(myJmol, `isosurface p4 hkl {${String(HInputs)} ${String(KInputs)} ${String(LInputs)}} boundbox;color isosurface red translucent`);
  HighlightBands(HInputs,KInputs,IInputs,LInputs);
}

init();
animate();
kinematic();

//console.log(SF)