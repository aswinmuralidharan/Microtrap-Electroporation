function out = model
%
% Microtrap_array_2D.m
%
% Model exported on Nov 29 2021, 18:35 by COMSOL 5.6.0.401.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('');

model.label('Array_2.mph');

model.param.set('Vapp', '90 [V]', 'Applied Voltage');
model.param.set('l', '1 [mm]', 'Distance between electrodes');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);

model.component('comp1').curvedInterior(false);

model.result.table.create('evl2', 'Table');

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').lengthUnit('mm');
model.component('comp1').geom('geom1').create('imp1', 'Import');
model.component('comp1').geom('geom1').feature('imp1').set('type', 'dxf');
model.component('comp1').geom('geom1').feature('imp1').set('filename', 'Trap_CAD_file.dxf');
model.component('comp1').geom('geom1').feature('imp1').set('alllayers', {'0'});
model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('r1').set('pos', {'-0.0270' '-0.35'});
model.component('comp1').geom('geom1').feature('r1').set('size', [1 1.4]);
model.component('comp1').geom('geom1').create('c1', 'Circle');
model.component('comp1').geom('geom1').feature('c1').active(false);
model.component('comp1').geom('geom1').feature('c1').set('pos', [0.4973 1.1]);
model.component('comp1').geom('geom1').feature('c1').set('r', 0.2);
model.component('comp1').geom('geom1').create('c2', 'Circle');
model.component('comp1').geom('geom1').feature('c2').active(false);
model.component('comp1').geom('geom1').feature('c2').set('pos', [0.4973 -0.4]);
model.component('comp1').geom('geom1').feature('c2').set('r', 0.2);
model.component('comp1').geom('geom1').create('sel1', 'ExplicitSelection');
model.component('comp1').geom('geom1').feature('sel1').label('Traps');
model.component('comp1').geom('geom1').feature('sel1').selection('selection').set('imp1(1)', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271 272 273 274 275 276 277 278 279 280 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300 301 302 303 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323 324 325 326 327 328 329 330 331 332 333 334 335 336 337 338 339 340 341 342 343 344 345 346 347 348 349 350 351 352 353 354 355 356 357 358 359 360 361 362 363 364 365 366 367 368 369 370 371 372 373 374 375 376 377 378 379 380 381 382 383 384 385 386 387 388 389 390 391 392 393 394 395 396 397 398 399 400 401 402 403 404 405 406]);
model.component('comp1').geom('geom1').create('dif1', 'Difference');
model.component('comp1').geom('geom1').feature('dif1').active(false);
model.component('comp1').geom('geom1').feature('dif1').selection('input').set({'r1'});
model.component('comp1').geom('geom1').feature('dif1').selection('input2').set({'imp1'});
model.component('comp1').geom('geom1').run;

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material.create('mat2', 'Common');
model.component('comp1').material('mat1').info.create('Composition');
model.component('comp1').material('mat1').propertyGroup('def').func.create('k_liquid_2', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('alpha_liquid_2', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('C_liquid_2', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('HC_liquid_2', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('VP_liquid_2', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('rho_liquid_2', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('TD_liquid_2', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('eta_liquid_1', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup.create('ThermalExpansion', 'Thermal expansion');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').func.create('dL_liquid_2', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').func.create('CTE_liquid_2', 'Piecewise');
model.component('comp1').material('mat2').selection.named('geom1_sel1');
model.component('comp1').material('mat2').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');

model.component('comp1').physics.create('ec', 'ConductiveMedia', 'geom1');
model.component('comp1').physics('ec').create('pot1', 'ElectricPotential', 1);
model.component('comp1').physics('ec').feature('pot1').selection.set([3]);
model.component('comp1').physics('ec').create('pot2', 'ElectricPotential', 1);
model.component('comp1').physics('ec').feature('pot2').selection.set([2]);

model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri2', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftri1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri1').selection.set([1]);
model.component('comp1').mesh('mesh1').feature('ftri1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri2').selection.named('geom1_sel1');
model.component('comp1').mesh('mesh1').feature('ftri2').create('size1', 'Size');

model.result.table('evl2').label('Evaluation 2D');
model.result.table('evl2').comments('Interactive 2D values');

model.component('comp1').view('view1').axis.set('xmin', 0.2929350733757019);
model.component('comp1').view('view1').axis.set('xmax', 0.6530649662017822);
model.component('comp1').view('view1').axis.set('ymin', 0.25232040882110596);
model.component('comp1').view('view1').axis.set('ymax', 0.44767963886260986);

model.component('comp1').material('mat1').label('H2O (water) [liquid]');
model.component('comp1').material('mat1').set('info', {'Composition' '' '11.2 H, 88.8 O (wt%)'});
model.component('comp1').material('mat1').propertyGroup('def').func('k_liquid_2').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('k_liquid_2').set('pieces', {'275.0' '370.0' '-0.9003748+0.008387698*T^1-1.118205E-5*T^2'});
model.component('comp1').material('mat1').propertyGroup('def').func('alpha_liquid_2').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('alpha_liquid_2').set('pieces', {'273.0' '283.0' '0.01032507-7.62815E-5*T^1+1.412474E-7*T^2'; '283.0' '373.0' '-0.002464185+1.947611E-5*T^1-5.049672E-8*T^2+4.616995E-11*T^3'});
model.component('comp1').material('mat1').propertyGroup('def').func('C_liquid_2').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('C_liquid_2').set('pieces', {'293.0' '373.0' '4035.841+0.492312*T^1'});
model.component('comp1').material('mat1').propertyGroup('def').func('HC_liquid_2').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('HC_liquid_2').set('pieces', {'293.0' '373.0' '72.64512+0.008861616*T^1'});
model.component('comp1').material('mat1').propertyGroup('def').func('VP_liquid_2').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('VP_liquid_2').set('pieces', {'280.0' '600.0' '(exp((-2.005122e+03/T-5.565700e-01*log10(T)+9.898790e+00-1.111690e+07/T^3)*log(10.0)))*1.333200e+02'});
model.component('comp1').material('mat1').propertyGroup('def').func('rho_liquid_2').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('rho_liquid_2').set('pieces', {'273.0' '283.0' '972.7584+0.2084*T^1-4.0E-4*T^2'; '283.0' '373.0' '345.28+5.749816*T^1-0.0157244*T^2+1.264375E-5*T^3'});
model.component('comp1').material('mat1').propertyGroup('def').func('TD_liquid_2').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('TD_liquid_2').set('pieces', {'273.0' '333.0' '8.04E-8+2.0E-10*T^1'});
model.component('comp1').material('mat1').propertyGroup('def').func('eta_liquid_1').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('eta_liquid_1').set('pieces', {'265.0' '293.0' '5.948859-0.08236196*T^1+4.287142E-4*T^2-9.938045E-7*T^3+8.65316E-10*T^4'; '293.0' '353.0' '0.410191-0.004753985*T^1+2.079795E-5*T^2-4.061698E-8*T^3+2.983925E-11*T^4'; '353.0' '423.0' '0.03625638-3.265463E-4*T^1+1.127139E-6*T^2-1.75363E-9*T^3+1.033976E-12*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', {'k_liquid_2(T[1/K])[W/(m*K)]' '0' '0' '0' 'k_liquid_2(T[1/K])[W/(m*K)]' '0' '0' '0' 'k_liquid_2(T[1/K])[W/(m*K)]'});
model.component('comp1').material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', {'(alpha_liquid_2(T[1/K])[1/K]+(Tempref-293[K])*if(abs(T-Tempref)>1e-3,(alpha_liquid_2(T[1/K])[1/K]-alpha_liquid_2(Tempref[1/K])[1/K])/(T-Tempref),d(alpha_liquid_2(T[1/K])[1/K],T)))/(1+alpha_liquid_2(Tempref[1/K])[1/K]*(Tempref-293[K]))' '0' '0' '0' '(alpha_liquid_2(T[1/K])[1/K]+(Tempref-293[K])*if(abs(T-Tempref)>1e-3,(alpha_liquid_2(T[1/K])[1/K]-alpha_liquid_2(Tempref[1/K])[1/K])/(T-Tempref),d(alpha_liquid_2(T[1/K])[1/K],T)))/(1+alpha_liquid_2(Tempref[1/K])[1/K]*(Tempref-293[K]))' '0' '0' '0' '(alpha_liquid_2(T[1/K])[1/K]+(Tempref-293[K])*if(abs(T-Tempref)>1e-3,(alpha_liquid_2(T[1/K])[1/K]-alpha_liquid_2(Tempref[1/K])[1/K])/(T-Tempref),d(alpha_liquid_2(T[1/K])[1/K],T)))/(1+alpha_liquid_2(Tempref[1/K])[1/K]*(Tempref-293[K]))'});
model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', 'C_liquid_2(T[1/K])[J/(kg*K)]');
model.component('comp1').material('mat1').propertyGroup('def').set('HC', 'HC_liquid_2(T[1/K])[J/(mol*K)]');
model.component('comp1').material('mat1').propertyGroup('def').set('VP', 'VP_liquid_2(T[1/K])[Pa]');
model.component('comp1').material('mat1').propertyGroup('def').set('density', 'rho_liquid_2(T[1/K])[kg/m^3]');
model.component('comp1').material('mat1').propertyGroup('def').set('TD', 'TD_liquid_2(T[1/K])[m^2/s]');
model.component('comp1').material('mat1').propertyGroup('def').set('dynamicviscosity', 'eta_liquid_1(T[1/K])[Pa*s]');
model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'1.5' '0' '0' '0' '1.5' '0' '0' '0' '1.5'});
model.component('comp1').material('mat1').propertyGroup('def').set('relpermittivity', {'80.1' '0' '0' '0' '80.1' '0' '0' '0' '80.1'});
model.component('comp1').material('mat1').propertyGroup('def').addInput('temperature');
model.component('comp1').material('mat1').propertyGroup('def').addInput('strainreferencetemperature');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').func('dL_liquid_2').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').func('dL_liquid_2').set('pieces', {'273.0' '283.0' '0.008486158-6.947021E-5*T^1+1.333373E-7*T^2'; '283.0' '373.0' '0.2324466-0.002030447*T^1+5.510259E-6*T^2-4.395999E-9*T^3'});
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').func('CTE_liquid_2').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').func('CTE_liquid_2').set('pieces', {'273.0' '283.0' '-6.947321E-5+2.666746E-7*T^1'; '283.0' '293.0' '-0.01363715+8.893977E-5*T^1-1.43925E-7*T^2'; '293.0' '373.0' '-0.00203045+1.102052E-5*T^1-1.3188E-8*T^2'});
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').set('alphatan', '');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').set('dL', '');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').set('alphatan', {'CTE_liquid_2(T[1/K])[1/K]' '0' '0' '0' 'CTE_liquid_2(T[1/K])[1/K]' '0' '0' '0' 'CTE_liquid_2(T[1/K])[1/K]'});
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').set('dL', {'(dL_liquid_2(T[1/K])-dL_liquid_2(Tempref[1/K]))/(1+dL_liquid_2(Tempref[1/K]))' '0' '0' '0' '(dL_liquid_2(T[1/K])-dL_liquid_2(Tempref[1/K]))/(1+dL_liquid_2(Tempref[1/K]))' '0' '0' '0' '(dL_liquid_2(T[1/K])-dL_liquid_2(Tempref[1/K]))/(1+dL_liquid_2(Tempref[1/K]))'});
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').addInput('temperature');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').addInput('strainreferencetemperature');
model.component('comp1').material('mat2').label('PDMS - Polydimethylsiloxane');
model.component('comp1').material('mat2').propertyGroup('def').set('thermalexpansioncoefficient', {'9e-4[1/K]' '0' '0' '0' '9e-4[1/K]' '0' '0' '0' '9e-4[1/K]'});
model.component('comp1').material('mat2').propertyGroup('def').set('heatcapacity', '1460[J/(kg*K)]');
model.component('comp1').material('mat2').propertyGroup('def').set('relpermittivity', {'2.75' '0' '0' '0' '2.75' '0' '0' '0' '2.75'});
model.component('comp1').material('mat2').propertyGroup('def').set('density', '970[kg/m^3]');
model.component('comp1').material('mat2').propertyGroup('def').set('thermalconductivity', {'0.16[W/(m*K)]' '0' '0' '0' '0.16[W/(m*K)]' '0' '0' '0' '0.16[W/(m*K)]'});
model.component('comp1').material('mat2').propertyGroup('def').set('electricconductivity', {'4*10^-13' '0' '0' '0' '4*10^-13' '0' '0' '0' '4*10^-13'});
model.component('comp1').material('mat2').propertyGroup('Enu').set('youngsmodulus', '750[kPa]');
model.component('comp1').material('mat2').propertyGroup('Enu').set('poissonsratio', '0.49');

model.component('comp1').physics('ec').prop('PortSweepSettings').set('PortParamName', 'PortName');
model.component('comp1').physics('ec').feature('cucn1').set('minput_temperature_src', 'userdef');
model.component('comp1').physics('ec').feature('cucn1').set('minput_numberdensity', 0);
model.component('comp1').physics('ec').feature('pot1').set('V0', 'Vapp');

model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hauto', 3);
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hauto', 4);
model.component('comp1').mesh('mesh1').run;

model.study.create('std1');
model.study('std1').create('param', 'Parametric');
model.study('std1').create('stat', 'Stationary');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('p1', 'Parametric');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature.remove('fcDef');

model.result.create('pg1', 'PlotGroup2D');
model.result('pg1').selection.geom('geom1', 2);
model.result('pg1').selection.set([1]);
model.result('pg1').create('surf1', 'Surface');
model.result('pg1').feature('surf1').set('expr', 'sqrt(ec.Ex*ec.Ex+ec.Ey*ec.Ey)/Vapp*l');

model.study('std1').feature('param').set('pname', {'Vapp'});
model.study('std1').feature('param').set('plistarr', {'90, 180, 270, 320'});
model.study('std1').feature('param').set('punit', {'V'});

model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').label('Compile Equations: Stationary');
model.sol('sol1').feature('v1').label('Dependent Variables 1.1');
model.sol('sol1').feature('v1').set('clistctrl', {'p1'});
model.sol('sol1').feature('v1').set('cname', {'Vapp'});
model.sol('sol1').feature('v1').set('clist', {'90[V] 180[V] 270[V] 320[V]'});
model.sol('sol1').feature('s1').label('Stationary Solver 1.1');
model.sol('sol1').feature('s1').set('probesel', 'none');
model.sol('sol1').feature('s1').feature('dDef').label('Direct 1');
model.sol('sol1').feature('s1').feature('aDef').label('Advanced 1');
model.sol('sol1').feature('s1').feature('p1').label('Parametric 1.1');
model.sol('sol1').feature('s1').feature('p1').set('control', 'param');
model.sol('sol1').feature('s1').feature('p1').set('pname', {'Vapp'});
model.sol('sol1').feature('s1').feature('p1').set('plistarr', {'90, 180, 270, 320'});
model.sol('sol1').feature('s1').feature('p1').set('punit', {'V'});
model.sol('sol1').feature('s1').feature('p1').set('pcontinuationmode', 'no');
model.sol('sol1').feature('s1').feature('fc1').label('Fully Coupled 1.1');
model.sol('sol1').runAll;

model.result('pg1').label('Electric field amplification');
model.result('pg1').set('looplevel', [1]);
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').feature('surf1').set('rangecoloractive', true);
model.result('pg1').feature('surf1').set('rangecolormax', 3);
model.result('pg1').feature('surf1').set('resolution', 'normal');

out = model;
