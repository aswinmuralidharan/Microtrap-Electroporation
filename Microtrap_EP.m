function out = model
%
% Microtrap_EP.m
%
% Model exported on Jan 28 2022, 10:25 by COMSOL 5.6.0.401.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('');

model.label('Microtrap_diff_2_output_solved_LR.mph');

model.param.set('sigma_e', '1.5 [S/m]', 'Extracellular conductivity');
model.param.set('epsilon_e', '80', 'extracellular permittivity');
model.param.set('sigma_i', '0.5[S/m]', 'Intracellular conductivity');
model.param.set('epsilon_i', '80', 'Intracellular permittivity');
model.param.set('epsilon_cm', '5', 'Cell membrane permittivity');
model.param.set('sigma_cm', '3e-7[S/m]', 'Cell membrane conductivity');
model.param.set('Vep', '258 [mV]', 'Threshold electroporation voltage');
model.param.set('alpha', '1e9[1/(m^2*s)]', 'Electroporation parameter');
model.param.set('N0', '1.5e5[cm^-2]', 'Equilibrium pore density when ITV = 0 V');
model.param.set('q', '2.46', 'Electroporation constant');
model.param.set('F', '96485[A*s/mol]', 'Faraday constant');
model.param.set('T', '25[degC]', 'Temperature');
model.param.set('R', '8.341[J/(K*mol)]', 'Universal gas constant');
model.param.set('w', '2.65', 'Energy barrier within pore');
model.param.set('n', '0.15', 'Relative entrance length of pores');
model.param.set('rp', '0.76e-9[m]', 'Pore radius');
model.param.set('dcm', '5[nm]', 'Membrane thickness');
model.param.set('Uapp', '3.595 [V]');
model.param.set('del', '0.7 [mm]');
model.param.set('E', 'Uapp/del');
model.param.set('tpulse', '5[ms]');
model.param.set('trise', '1e-06[s]');
model.param.set('dmesh1', '5e-4');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 3);

model.component('comp1').curvedInterior(false);

model.result.table.create('evl3', 'Table');
model.result.table.create('tbl1', 'Table');
model.result.table.create('tbl2', 'Table');
model.result.table.create('tbl3', 'Table');
model.result.table.create('tbl4', 'Table');

model.func.create('an1', 'Analytic');
model.func('an1').set('funcname', 'pulse');
model.func('an1').set('expr', 'flc1hs(t-trise/2,trise/2) - flc1hs(t-tpulse-trise/2,trise/2)');
model.func('an1').set('args', {'t'});
model.func('an1').set('argunit', 's');
model.func('an1').set('plotargs', {'t' '0' '1'});

model.component('comp1').mesh.create('mesh1');
model.component('comp1').mesh.create('mesh2');

model.component('comp1').geom('geom1').lengthUnit('mm');
model.component('comp1').geom('geom1').create('sph1', 'Sphere');
model.component('comp1').geom('geom1').feature('sph1').set('r', '10 [um]');
model.component('comp1').geom('geom1').create('blk1', 'Block');
model.component('comp1').geom('geom1').feature('blk1').set('pos', {'0 [um]' '0' '0'});
model.component('comp1').geom('geom1').feature('blk1').set('base', 'center');
model.component('comp1').geom('geom1').feature('blk1').set('size', {'48.8 [um]' '48.8 [um]' '35 [um]'});
model.component('comp1').geom('geom1').create('blk2', 'Block');
model.component('comp1').geom('geom1').feature('blk2').set('pos', {'-50 [um]' '0[um]' '0 [um]'});
model.component('comp1').geom('geom1').feature('blk2').set('base', 'center');
model.component('comp1').geom('geom1').feature('blk2').set('size', {'100 [um]' '100 [um]' '100 [um]'});
model.component('comp1').geom('geom1').create('blk5', 'Block');
model.component('comp1').geom('geom1').feature('blk5').set('pos', {'0[um]' '0[um]' '-50 [um]'});
model.component('comp1').geom('geom1').feature('blk5').set('base', 'center');
model.component('comp1').geom('geom1').feature('blk5').set('size', {'100 [um]' '100 [um]' '100 [um]'});
model.component('comp1').geom('geom1').create('imp1', 'Import');
model.component('comp1').geom('geom1').feature('imp1').set('type', 'cad');
model.component('comp1').geom('geom1').feature('imp1').set('filename', 'single_trap_3.dxf');
model.component('comp1').geom('geom1').feature('imp1').set('importtol', 1.0E-7);
model.component('comp1').geom('geom1').feature('imp1').importData;
model.component('comp1').geom('geom1').create('mov1', 'Move');
model.component('comp1').geom('geom1').feature('mov1').set('displx', '0.1[um]');
model.component('comp1').geom('geom1').feature('mov1').set('disply', '0.0');
model.component('comp1').geom('geom1').feature('mov1').selection('input').set({'imp1'});
model.component('comp1').geom('geom1').create('blk3', 'Block');
model.component('comp1').geom('geom1').feature('blk3').set('pos', {'0' '0' '25 [um]'});
model.component('comp1').geom('geom1').feature('blk3').set('base', 'center');
model.component('comp1').geom('geom1').feature('blk3').set('size', {'100 [um]' '100 [um]' '15 [um]'});
model.component('comp1').geom('geom1').create('co1', 'Compose');
model.component('comp1').geom('geom1').feature('co1').set('formula', 'blk1+sph1+mov1-(blk2)-blk5-blk3');
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('blk1');

model.component('comp1').variable.create('var1');
model.component('comp1').variable('var1').set('ITV', 'V-V2+0.005', 'Induced Transmembrane Voltage');
model.component('comp1').variable('var1').set('sigma_p', '(sigma_e-sigma_i)/log(sigma_e/sigma_i)', 'Pore conductance');
model.component('comp1').variable.create('var2');
model.component('comp1').variable('var2').set('A', 'alpha*exp((ITV/Vep)^2)');
model.component('comp1').variable('var2').set('B', 'N0*exp(q*(ITV/Vep)^2)');
model.component('comp1').variable('var2').set('vm', 'ITV*F/(R*T)');
model.component('comp1').variable('var2').set('K0', 'w/((w*(1-2*n)+2*n)*exp(w)-2*n)');
model.component('comp1').variable('var2').set('K1', '(w*exp(w-n*vm)-n*vm)/(w-n*vm)');
model.component('comp1').variable('var2').set('K2', '(w*exp(w+n*vm)+n*vm)/(w+n*vm)');
model.component('comp1').variable('var2').set('inter', 'K0*(ITV==0)+(-1+exp(vm))/(K1*exp(vm)-K2+(ITV==0))');
model.component('comp1').variable('var2').set('sigma_ep', 'pi*rp^2*sigma_p*u*inter');
model.component('comp1').variable('var2').set('sigma_ep_new', 'u*(2*pi*rp^2*sigma_p*dcm/(pi*rp + 2*dcm))');
model.component('comp1').variable.create('var3');
model.component('comp1').variable('var3').set('ITVav', 'aveop1(ITV)');

model.component('comp1').view('view1').hideEntities.create('hide1');
model.component('comp1').view('view1').hideEntities('hide1').geom('geom1', 2);
model.component('comp1').view('view1').hideEntities('hide1').set([10 11 12 13 14 15 16 17 18 19 20 21 22 23]);

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material.create('mat2', 'Common');
model.component('comp1').material.create('mat3', 'Common');
model.component('comp1').material.create('mat4', 'Common');
model.component('comp1').material('mat1').selection.set([1]);
model.component('comp1').material('mat2').selection.set([2]);
model.component('comp1').material('mat3').selection.geom('geom1', 2);
model.component('comp1').material('mat3').selection.set([6 8]);
model.component('comp1').material('mat4').selection.set([3]);
model.component('comp1').material('mat4').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');

model.component('comp1').cpl.create('aveop1', 'Average');
model.component('comp1').cpl('aveop1').selection.geom('geom1', 1);
model.component('comp1').cpl('aveop1').selection.set([6 9]);

model.component('comp1').physics.create('ec', 'ConductiveMedia', 'geom1');
model.component('comp1').physics('ec').selection.set([1 3]);
model.component('comp1').physics('ec').create('pot1', 'ElectricPotential', 2);
model.component('comp1').physics('ec').feature('pot1').selection.set([9]);
model.component('comp1').physics('ec').create('dimp1', 'DistributedImpedance', 2);
model.component('comp1').physics('ec').feature('dimp1').selection.set([6 8]);
model.component('comp1').physics('ec').create('gnd1', 'Ground', 2);
model.component('comp1').physics('ec').feature('gnd1').selection.set([2]);
model.component('comp1').physics.create('ec2', 'ConductiveMedia', 'geom1');
model.component('comp1').physics('ec2').selection.set([2]);
model.component('comp1').physics('ec2').create('dimp1', 'DistributedImpedance', 2);
model.component('comp1').physics('ec2').feature('dimp1').selection.set([6 8]);
model.component('comp1').physics.create('wb', 'WeakFormBoundaryPDE', 'geom1');
model.component('comp1').physics('wb').prop('Units').set('DependentVariableQuantity', 'none');
model.component('comp1').physics('wb').prop('Units').set('CustomDependentVariableUnit', 'm^-2');
model.component('comp1').physics('wb').selection.set([6 8]);

model.component('comp1').mesh('mesh1').create('ftet2', 'FreeTet');
model.component('comp1').mesh('mesh1').feature('ftet2').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').selection.set([6 8]);

model.result.table('evl3').label('Evaluation 3D');
model.result.table('evl3').comments('Interactive 3D values');
model.result.table('tbl3').label('Table 3');

model.component('comp1').view('view1').set('showgrid', false);
model.component('comp1').view('view1').set('scenelight', false);
model.component('comp1').view('view1').set('showselection', false);
model.component('comp1').view('view1').set('showmaterial', true);

model.component('comp1').material('mat1').set('color', 'gray');
model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'sigma_e' '0' '0' '0' 'sigma_e' '0' '0' '0' 'sigma_e'});
model.component('comp1').material('mat1').propertyGroup('def').set('relpermittivity', {'epsilon_e' '0' '0' '0' 'epsilon_e' '0' '0' '0' 'epsilon_e'});
model.component('comp1').material('mat2').set('color', 'green');
model.component('comp1').material('mat2').propertyGroup('def').set('electricconductivity', {'sigma_i' '0' '0' '0' 'sigma_i' '0' '0' '0' 'sigma_i'});
model.component('comp1').material('mat2').propertyGroup('def').set('relpermittivity', {'epsilon_i' '0' '0' '0' 'epsilon_i' '0' '0' '0' 'epsilon_i'});
model.component('comp1').material('mat3').set('color', 'black');
model.component('comp1').material('mat3').propertyGroup('def').set('electricconductivity', '');
model.component('comp1').material('mat3').propertyGroup('def').set('relpermittivity', '');
model.component('comp1').material('mat3').propertyGroup('def').set('electricconductivity', {'sigma_cm+sigma_ep' '0' '0' '0' 'sigma_cm+sigma_ep' '0' '0' '0' 'sigma_cm+sigma_ep'});
model.component('comp1').material('mat3').propertyGroup('def').set('relpermittivity', {'epsilon_cm' '0' '0' '0' 'epsilon_cm' '0' '0' '0' 'epsilon_cm'});
model.component('comp1').material('mat4').label('PDMS - Polydimethylsiloxane');
model.component('comp1').material('mat4').propertyGroup('def').set('thermalexpansioncoefficient', {'9e-4[1/K]' '0' '0' '0' '9e-4[1/K]' '0' '0' '0' '9e-4[1/K]'});
model.component('comp1').material('mat4').propertyGroup('def').set('heatcapacity', '1460[J/(kg*K)]');
model.component('comp1').material('mat4').propertyGroup('def').set('relpermittivity', {'2.75' '0' '0' '0' '2.75' '0' '0' '0' '2.75'});
model.component('comp1').material('mat4').propertyGroup('def').set('density', '970[kg/m^3]');
model.component('comp1').material('mat4').propertyGroup('def').set('thermalconductivity', {'0.16[W/(m*K)]' '0' '0' '0' '0.16[W/(m*K)]' '0' '0' '0' '0.16[W/(m*K)]'});
model.component('comp1').material('mat4').propertyGroup('def').set('electricconductivity', {'4e-13' '0' '0' '0' '4e-13' '0' '0' '0' '4e-13'});
model.component('comp1').material('mat4').propertyGroup('Enu').set('youngsmodulus', '750[kPa]');
model.component('comp1').material('mat4').propertyGroup('Enu').set('poissonsratio', '0.49');

model.component('comp1').physics('ec').prop('PortSweepSettings').set('PortParamName', 'PortName');
model.component('comp1').physics('ec').feature('cucn1').set('minput_temperature_src', 'userdef');
model.component('comp1').physics('ec').feature('cucn1').set('minput_numberdensity', 0);
model.component('comp1').physics('ec').feature('pot1').set('V0', 'Uapp*pulse(t)');
model.component('comp1').physics('ec').feature('dimp1').set('Vref', 'V2');
model.component('comp1').physics('ec').feature('dimp1').set('ds', 'dcm');
model.component('comp1').physics('ec2').prop('PortSweepSettings').set('PortParamName', 'PortName');
model.component('comp1').physics('ec2').feature('cucn1').set('minput_temperature_src', 'userdef');
model.component('comp1').physics('ec2').feature('cucn1').set('minput_numberdensity', 0);
model.component('comp1').physics('ec2').feature('dimp1').set('Vref', 'V');
model.component('comp1').physics('ec2').feature('dimp1').set('ds', 'dcm');
model.component('comp1').physics('wb').prop('Units').set('CustomSourceTermUnit', 'm^-2*s^-1');
model.component('comp1').physics('wb').feature('wfeq1').set('weak', 'test(u)*(A-A/B*u-ut)');

model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hmax', 'dmesh1');
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').run;

model.study.create('std1');
model.study('std1').create('param', 'Parametric');
model.study('std1').create('time', 'Transient');
model.study.create('std2');
model.study('std2').create('param', 'Parametric');
model.study('std2').create('time', 'Transient');
model.study.create('std3');
model.study('std3').create('param', 'Parametric');
model.study('std3').create('time', 'Transient');

model.sol.create('sol6');
model.sol('sol6').study('std2');
model.sol('sol6').attach('std2');
model.sol('sol6').create('st1', 'StudyStep');
model.sol('sol6').create('v1', 'Variables');
model.sol('sol6').create('t1', 'Time');
model.sol('sol6').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol6').feature('t1').feature.remove('fcDef');
model.sol.create('sol7');
model.sol('sol7').study('std2');
model.sol('sol7').label('Parametric Solutions 2');
model.sol.create('sol12');
model.sol('sol12').study('std1');
model.sol('sol12').attach('std1');
model.sol('sol12').create('st1', 'StudyStep');
model.sol('sol12').create('v1', 'Variables');
model.sol('sol12').create('t1', 'Time');
model.sol('sol12').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol12').feature('t1').feature.remove('fcDef');
model.sol.create('sol13');
model.sol('sol13').study('std1');
model.sol('sol13').label('Parametric Solutions 1');
model.sol.create('sol17');
model.sol('sol17').study('std3');
model.sol('sol17').attach('std3');
model.sol('sol17').create('st1', 'StudyStep');
model.sol('sol17').create('v1', 'Variables');
model.sol('sol17').create('t1', 'Time');
model.sol('sol17').feature('t1').create('taDef', 'TimeAdaption');
model.sol('sol17').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol17').feature('t1').feature.remove('fcDef');
model.sol.create('sol18');
model.sol('sol18').study('std3');
model.sol.create('sol19');
model.sol('sol19').study('std3');
model.sol('sol19').label('Parametric Solutions 3');

model.batch.create('p2', 'Parametric');
model.batch.create('p3', 'Parametric');
model.batch.create('p4', 'Parametric');
model.batch('p2').create('so1', 'Solutionseq');
model.batch('p3').create('so1', 'Solutionseq');
model.batch('p4').create('so1', 'Solutionseq');
model.batch('p2').study('std2');
model.batch('p3').study('std1');
model.batch('p4').study('std3');

model.result.dataset.create('dset8', 'Solution');
model.result.dataset.create('dset9', 'Solution');
model.result.dataset('dset3').set('solution', 'sol6');
model.result.dataset('dset4').set('solution', 'sol7');
model.result.dataset('dset5').set('solution', 'sol12');
model.result.dataset('dset6').set('solution', 'sol13');
model.result.dataset('dset7').set('solution', 'sol17');
model.result.dataset('dset8').set('solution', 'sol18');
model.result.dataset('dset9').set('solution', 'sol19');
model.result.dataset.remove('dset1');
model.result.dataset.remove('dset2');
model.result.create('pg18', 'PlotGroup1D');
model.result.create('pg20', 'PlotGroup1D');
model.result.create('pg14', 'PlotGroup1D');
model.result.create('pg19', 'PlotGroup1D');
model.result.create('pg21', 'PlotGroup3D');
model.result('pg18').set('data', 'dset6');
model.result('pg18').create('lngr1', 'LineGraph');
model.result('pg18').create('lngr2', 'LineGraph');
model.result('pg18').feature('lngr1').selection.set([6 8 9 12]);
model.result('pg18').feature('lngr1').set('expr', 'u');
model.result('pg18').feature('lngr2').set('data', 'dset9');
model.result('pg18').feature('lngr2').selection.set([6 8 9 12]);
model.result('pg18').feature('lngr2').set('expr', 'u');
model.result('pg20').set('data', 'dset9');
model.result('pg20').create('lngr1', 'LineGraph');
model.result('pg20').feature('lngr1').selection.set([6 8 9 12]);
model.result('pg20').feature('lngr1').set('expr', 'u');
model.result('pg14').set('data', 'dset4');
model.result('pg14').create('lngr1', 'LineGraph');
model.result('pg14').feature('lngr1').selection.set([6 8 9 12]);
model.result('pg14').feature('lngr1').set('expr', 'u');
model.result('pg19').create('ptgr1', 'PointGraph');
model.result('pg19').feature('ptgr1').selection.set([3 5]);
model.result('pg19').feature('ptgr1').set('expr', 'abs(ITV)');
model.result('pg21').set('data', 'dset6');
model.result('pg21').create('surf1', 'Surface');
model.result('pg21').feature('surf1').set('expr', 'ITV');
model.result.report.create('rpt1', 'Report');
model.result.report('rpt1').create('tp1', 'TitlePage');
model.result.report('rpt1').create('toc1', 'TableOfContents');
model.result.report('rpt1').create('sec1', 'Section');
model.result.report('rpt1').create('sec2', 'Section');
model.result.report('rpt1').create('sec3', 'Section');
model.result.report('rpt1').create('sec4', 'Section');
model.result.report('rpt1').feature('sec1').create('root1', 'Model');
model.result.report('rpt1').feature('sec1').create('sec1', 'Section');
model.result.report('rpt1').feature('sec1').create('sec2', 'Section');
model.result.report('rpt1').feature('sec1').feature('sec1').create('param1', 'Parameter');
model.result.report('rpt1').feature('sec1').feature('sec2').create('sec1', 'Section');
model.result.report('rpt1').feature('sec1').feature('sec2').feature('sec1').create('func1', 'Functions');
model.result.report('rpt1').feature('sec2').create('comp1', 'ModelNode');
model.result.report('rpt1').feature('sec2').create('sec1', 'Section');
model.result.report('rpt1').feature('sec2').create('sec2', 'Section');
model.result.report('rpt1').feature('sec2').create('sec3', 'Section');
model.result.report('rpt1').feature('sec2').create('sec4', 'Section');
model.result.report('rpt1').feature('sec2').create('sec5', 'Section');
model.result.report('rpt1').feature('sec2').create('sec6', 'Section');
model.result.report('rpt1').feature('sec2').create('sec7', 'Section');
model.result.report('rpt1').feature('sec2').feature('sec1').create('sec1', 'Section');
model.result.report('rpt1').feature('sec2').feature('sec1').create('sec2', 'Section');
model.result.report('rpt1').feature('sec2').feature('sec1').create('sec3', 'Section');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec1').create('sec1', 'Section');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec1').create('sec2', 'Section');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec1').create('sec3', 'Section');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec1').feature('sec1').create('var1', 'Variables');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec1').feature('sec2').create('var1', 'Variables');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec1').feature('sec3').create('var1', 'Variables');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec2').create('sec1', 'Section');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec2').feature('sec1').create('cpl1', 'ComponentCoupling');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec3').create('sec1', 'Section');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec3').feature('sec1').create('csys1', 'CoordinateSystem');
model.result.report('rpt1').feature('sec2').feature('sec2').create('geom1', 'Geometry');
model.result.report('rpt1').feature('sec2').feature('sec3').create('sec1', 'Section');
model.result.report('rpt1').feature('sec2').feature('sec3').create('sec2', 'Section');
model.result.report('rpt1').feature('sec2').feature('sec3').create('sec3', 'Section');
model.result.report('rpt1').feature('sec2').feature('sec3').feature('sec1').create('mat1', 'Material');
model.result.report('rpt1').feature('sec2').feature('sec3').feature('sec2').create('mat1', 'Material');
model.result.report('rpt1').feature('sec2').feature('sec3').feature('sec3').create('mat1', 'Material');
model.result.report('rpt1').feature('sec2').feature('sec4').create('phys1', 'Physics');
model.result.report('rpt1').feature('sec2').feature('sec5').create('phys1', 'Physics');
model.result.report('rpt1').feature('sec2').feature('sec6').create('phys1', 'Physics');
model.result.report('rpt1').feature('sec2').feature('sec7').create('mesh1', 'Mesh');
model.result.report('rpt1').feature('sec3').create('std1', 'Study');
model.result.report('rpt1').feature('sec4').create('sec1', 'Section');
model.result.report('rpt1').feature('sec4').create('sec2', 'Section');
model.result.report('rpt1').feature('sec4').create('sec3', 'Section');
model.result.report('rpt1').feature('sec4').feature('sec1').create('sec1', 'Section');
model.result.report('rpt1').feature('sec4').feature('sec1').create('sec2', 'Section');
model.result.report('rpt1').feature('sec4').feature('sec1').feature('sec1').create('dset1', 'DataSet');
model.result.report('rpt1').feature('sec4').feature('sec1').feature('sec2').create('dset1', 'DataSet');
model.result.report('rpt1').feature('sec4').feature('sec2').create('sec1', 'Section');
model.result.report('rpt1').feature('sec4').feature('sec2').feature('sec1').create('mtbl1', 'Table');
model.result.report('rpt1').feature('sec4').feature('sec3').create('sec1', 'Section');
model.result.report('rpt1').feature('sec4').feature('sec3').create('sec2', 'Section');
model.result.report('rpt1').feature('sec4').feature('sec3').create('sec3', 'Section');
model.result.report('rpt1').feature('sec4').feature('sec3').create('sec4', 'Section');
model.result.report('rpt1').feature('sec4').feature('sec3').create('sec5', 'Section');
model.result.report('rpt1').feature('sec4').feature('sec3').feature('sec1').create('pg1', 'PlotGroup');
model.result.report('rpt1').feature('sec4').feature('sec3').feature('sec2').create('pg1', 'PlotGroup');
model.result.report('rpt1').feature('sec4').feature('sec3').feature('sec3').create('pg1', 'PlotGroup');
model.result.report('rpt1').feature('sec4').feature('sec3').feature('sec4').create('pg1', 'PlotGroup');
model.result.report('rpt1').feature('sec4').feature('sec3').feature('sec5').create('pg1', 'PlotGroup');

model.study('std1').feature('param').set('pname', {'Uapp'});
model.study('std1').feature('param').set('plistarr', {'3.595, 7.3644, 11.047'});
model.study('std1').feature('param').set('punit', {'V'});
model.study('std1').feature('time').set('tlist', '10^{range(log10(1.0e-8),(log10(tpulse[1/s])+8)/100,log10(tpulse[1/s]))}');
model.study('std1').feature('time').set('usertol', true);
model.study('std1').feature('time').set('rtol', '0.0001');
model.study('std1').feature('time').set('mesh', {'geom1' 'mesh1'});
model.study('std2').label('Study 2 MeshSize');
model.study('std2').feature('param').set('pname', {'dmesh1'});
model.study('std2').feature('param').set('plistarr', {'1e-3,5e-4,2.5e-4,1e-4'});
model.study('std2').feature('param').set('punit', {''});
model.study('std2').feature('time').set('tlist', '10^{range(log10(1e-8),(log10(tpulse[1/s])+8)/100,log10(tpulse[1/s]))}');
model.study('std2').feature('time').set('usertol', true);
model.study('std2').feature('time').set('rtol', '0.0001');
model.study('std2').feature('time').set('mesh', {'geom1' 'mesh1'});
model.study('std3').label('Study 3 AdaptiveMesh');
model.study('std3').feature('param').set('pname', {'Uapp'});
model.study('std3').feature('param').set('plistarr', {'3.595, 7.3644, 11.047'});
model.study('std3').feature('param').set('punit', {'V'});
model.study('std3').feature('time').set('tlist', '10^{range(log10(1.0e-8),(log10(tpulse[1/s])+8)/100,log10(tpulse[1/s]))}');
model.study('std3').feature('time').set('usertol', true);
model.study('std3').feature('time').set('rtol', '0.0001');
model.study('std3').feature('time').set('timeadaption', true);

model.sol('sol6').attach('std2');
model.sol('sol6').feature('st1').label('Compile Equations: Time Dependent');
model.sol('sol6').feature('v1').label('Dependent Variables 1.1');
model.sol('sol6').feature('v1').set('clist', {'10^{range(log10(1e-8),(log10(tpulse[1/s])+8)/100,log10(tpulse[1/s]))}' '1e-11[s]'});
model.sol('sol6').feature('t1').label('Time-Dependent Solver 1.1');
model.sol('sol6').feature('t1').set('control', 'user');
model.sol('sol6').feature('t1').set('tlist', '10^{range(log10(1e-8),(log10(tpulse[1/s])+8)/100,log10(tpulse[1/s]))}');
model.sol('sol6').feature('t1').set('rtol', '0.0001');
model.sol('sol6').feature('t1').set('tstepsbdf', 'strict');
model.sol('sol6').feature('t1').set('initialstepbdf', '1e-11');
model.sol('sol6').feature('t1').set('initialstepbdfactive', true);
model.sol('sol6').feature('t1').feature('dDef').label('Direct 1');
model.sol('sol6').feature('t1').feature('aDef').label('Advanced 1');
model.sol('sol6').feature('t1').feature('fc1').label('Fully Coupled 1.1');
model.sol('sol6').runAll;
model.sol('sol12').attach('std1');
model.sol('sol12').feature('st1').label('Compile Equations: Time Dependent');
model.sol('sol12').feature('v1').label('Dependent Variables 1.1');
model.sol('sol12').feature('v1').set('clist', {'10^{range(log10(1.0e-8),(log10(tpulse[1/s])+8)/100,log10(tpulse[1/s]))}' '1e-9[s]'});
model.sol('sol12').feature('t1').label('Time-Dependent Solver 1.1');
model.sol('sol12').feature('t1').set('tlist', '10^{range(log10(1.0e-8),(log10(tpulse[1/s])+8)/100,log10(tpulse[1/s]))}');
model.sol('sol12').feature('t1').set('rtol', '0.0001');
model.sol('sol12').feature('t1').set('tstepsbdf', 'strict');
model.sol('sol12').feature('t1').set('initialstepbdf', '1e-9');
model.sol('sol12').feature('t1').set('initialstepbdfactive', true);
model.sol('sol12').feature('t1').feature('dDef').label('Direct 1');
model.sol('sol12').feature('t1').feature('aDef').label('Advanced 1');
model.sol('sol12').feature('t1').feature('fc1').label('Fully Coupled 1.1');
model.sol('sol12').runAll;
model.sol('sol17').attach('std3');
model.sol('sol17').feature('st1').label('Compile Equations: Time Dependent');
model.sol('sol17').feature('v1').label('Dependent Variables 1.1');
model.sol('sol17').feature('v1').set('clist', {'10^{range(log10(1.0e-8),(log10(tpulse[1/s])+8)/100,log10(tpulse[1/s]))}' '1e-9[s]'});
model.sol('sol17').feature('t1').label('Time-Dependent Solver 1.1');
model.sol('sol17').feature('t1').set('tlist', '10^{range(log10(1.0e-8),(log10(tpulse[1/s])+8)/100,log10(tpulse[1/s]))}');
model.sol('sol17').feature('t1').set('rtol', '0.0001');
model.sol('sol17').feature('t1').set('tstepsbdf', 'strict');
model.sol('sol17').feature('t1').set('initialstepbdf', '1e-9');
model.sol('sol17').feature('t1').set('initialstepbdfactive', true);
model.sol('sol17').feature('t1').feature('taDef').label('Adaptive Mesh Refinement 1');
model.sol('sol17').feature('t1').feature('taDef').set('tadapsol', 'sol18');
model.sol('sol17').feature('t1').feature('taDef').set('tadapmesh', {'mesh3' 'mesh4' 'mesh5' 'mesh6' 'mesh7' 'mesh8' 'mesh9' 'mesh10' 'mesh11' 'mesh12'});
model.sol('sol17').feature('t1').feature('dDef').label('Direct 1');
model.sol('sol17').feature('t1').feature('aDef').label('Advanced 1');
model.sol('sol17').feature('t1').feature('fc1').label('Fully Coupled 1.1');
model.sol('sol17').runAll;
model.sol('sol18').label('Refined Mesh Solution 1');

model.batch('p2').set('pname', {'dmesh1'});
model.batch('p2').set('plistarr', {'1e-3,5e-4,2.5e-4,1e-4'});
model.batch('p2').set('punit', {''});
model.batch('p2').set('err', true);
model.batch('p2').feature('so1').set('seq', 'sol6');
model.batch('p2').feature('so1').set('store', true);
model.batch('p2').feature('so1').set('psol', 'sol7');
model.batch('p2').feature('so1').set('param', {'"dmesh1","0.001"' '"dmesh1","5E-4"' '"dmesh1","2.5E-4"' '"dmesh1","1E-4"'});
model.batch('p2').attach('std2');
model.batch('p2').run;
model.batch('p3').set('control', 'param');
model.batch('p3').set('pname', {'Uapp'});
model.batch('p3').set('plistarr', {'3.595, 7.3644, 11.047'});
model.batch('p3').set('punit', {'V'});
model.batch('p3').set('err', true);
model.batch('p3').feature('so1').set('seq', 'sol12');
model.batch('p3').feature('so1').set('psol', 'sol13');
model.batch('p3').feature('so1').set('param', {'"Uapp","3.595"' '"Uapp","7.3644"' '"Uapp","11.047"'});
model.batch('p3').attach('std1');
model.batch('p3').run;
model.batch('p4').set('control', 'param');
model.batch('p4').set('pname', {'Uapp'});
model.batch('p4').set('plistarr', {'3.595, 7.3644, 11.047'});
model.batch('p4').set('punit', {'V'});
model.batch('p4').set('err', true);
model.batch('p4').feature('so1').set('seq', 'sol17');
model.batch('p4').feature('so1').set('psol', 'sol19');
model.batch('p4').feature('so1').set('param', {'"Uapp","3.595"' '"Uapp","7.3644"' '"Uapp","11.047"'});
model.batch('p4').attach('std3');
model.batch('p4').run;

model.result('pg18').label('Result-Constantmesh');
model.result('pg18').set('looplevelinput', {'last' 'all'});
model.result('pg18').set('xlabel', 'Arc length (mm)');
model.result('pg18').set('ylabel', 'Dependent variable u (1/m<sup>2</sup>)');
model.result('pg18').set('xlabelactive', false);
model.result('pg18').set('ylabelactive', false);
model.result('pg18').feature('lngr1').set('linewidth', 4);
model.result('pg18').feature('lngr1').set('legend', true);
model.result('pg18').feature('lngr1').set('legendmethod', 'manual');
model.result('pg18').feature('lngr1').set('legends', {'3.595 V' '7.3644 V' '11.047 V'});
model.result('pg18').feature('lngr1').set('resolution', 'normal');
model.result('pg18').feature('lngr2').label('AdaptiveMesh');
model.result('pg18').feature('lngr2').set('looplevelinput', {'last' 'all'});
model.result('pg18').feature('lngr2').set('linestyle', 'dashed');
model.result('pg18').feature('lngr2').set('linecolor', 'black');
model.result('pg18').feature('lngr2').set('linewidth', 2);
model.result('pg18').feature('lngr2').set('legend', true);
model.result('pg18').feature('lngr2').set('legendmethod', 'manual');
model.result('pg18').feature('lngr2').set('legends', {'3.595 V' '7.3644 V' '11.047 V'});
model.result('pg18').feature('lngr2').set('resolution', 'normal');
model.result('pg20').label('Result-AdaptiveMesh ');
model.result('pg20').set('looplevelinput', {'last' 'all'});
model.result('pg20').set('xlabel', 'Arc length (mm)');
model.result('pg20').set('ylabel', 'Dependent variable u (1/m<sup>2</sup>)');
model.result('pg20').set('xlabelactive', false);
model.result('pg20').set('ylabelactive', false);
model.result('pg20').feature('lngr1').set('legend', true);
model.result('pg20').feature('lngr1').set('legendmethod', 'manual');
model.result('pg20').feature('lngr1').set('legends', {'3.595 V' '7.3644 V' '11.047 V'});
model.result('pg20').feature('lngr1').set('resolution', 'normal');
model.result('pg14').label('Influence-MeshSize');
model.result('pg14').set('looplevelinput', {'last' 'all'});
model.result('pg14').set('xlabel', 'Arc length (mm)');
model.result('pg14').set('ylabel', 'Dependent variable u (1/m<sup>2</sup>)');
model.result('pg14').set('xlabelactive', false);
model.result('pg14').set('ylabelactive', false);
model.result('pg14').feature('lngr1').set('legend', true);
model.result('pg14').feature('lngr1').set('legendmethod', 'manual');
model.result('pg14').feature('lngr1').set('legends', {'dmesh1=1.0 [um]' 'dmesh1= 0.5 [um]' 'dmesh1= 0.25 [um]' 'dmesh1=0.1 [um]'});
model.result('pg14').feature('lngr1').set('resolution', 'normal');
model.result('pg19').label('ITV');
model.result('pg19').set('xlabel', 'Time (s)');
model.result('pg19').set('ylabel', 'abs(ITV) (V)');
model.result('pg19').set('xlabelactive', false);
model.result('pg19').set('ylabelactive', false);
model.result('pg21').label('ITV SurfPlot');
model.result('pg21').feature('surf1').set('resolution', 'normal');
model.result.report('rpt1').set('templatesource', 'brief');
model.result.report('rpt1').set('format', 'docx');
model.result.report('rpt1').set('filename', 'Untitled.docx');
model.result.report('rpt1').feature('tp1').label('Microtrap electroporation');
model.result.report('rpt1').feature('tp1').set('frontmatterlayout', 'headings');
model.result.report('rpt1').feature('sec1').label('Global Definitions');
model.result.report('rpt1').feature('sec1').feature('sec1').label('Parameters');
model.result.report('rpt1').feature('sec1').feature('sec2').label('Functions');
model.result.report('rpt1').feature('sec1').feature('sec2').feature('sec1').label('Analytic 1');
model.result.report('rpt1').feature('sec1').feature('sec2').feature('sec1').set('source', 'firstchild');
model.result.report('rpt1').feature('sec1').feature('sec2').feature('sec1').feature('func1').set('includesettings', false);
model.result.report('rpt1').feature('sec2').label('Component 1');
model.result.report('rpt1').feature('sec2').feature('comp1').set('includeunitsystem', false);
model.result.report('rpt1').feature('sec2').feature('sec1').label('Definitions');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec1').label('Variables');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec1').feature('sec1').label('Variables 1');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec1').feature('sec1').set('source', 'firstchild');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec1').feature('sec2').label('Variables 2');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec1').feature('sec2').set('source', 'firstchild');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec1').feature('sec2').feature('var1').set('noderef', 'var2');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec1').feature('sec2').feature('var1').set('children', {'A' 'on';  ...
'B' 'on';  ...
'vm' 'on';  ...
'K0' 'on';  ...
'K1' 'on';  ...
'K2' 'on';  ...
'inter' 'on';  ...
'sigma_ep' 'on';  ...
'sigma_ep_new' 'on'});
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec1').feature('sec3').label('Variables 3');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec1').feature('sec3').set('source', 'firstchild');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec1').feature('sec3').feature('var1').set('noderef', 'var3');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec1').feature('sec3').feature('var1').set('children', {'ITVav' 'on'});
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec2').label('Nonlocal Couplings');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec2').feature('sec1').label('Average 1');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec2').feature('sec1').set('source', 'firstchild');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec2').feature('sec1').feature('cpl1').set('includesettings', false);
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec3').label('Coordinate Systems');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec3').feature('sec1').label('Boundary System 1');
model.result.report('rpt1').feature('sec2').feature('sec1').feature('sec3').feature('sec1').set('source', 'firstchild');
model.result.report('rpt1').feature('sec2').feature('sec2').label('Geometry 1');
model.result.report('rpt1').feature('sec2').feature('sec2').set('source', 'firstchild');
model.result.report('rpt1').feature('sec2').feature('sec2').feature('geom1').set('includestats', false);
model.result.report('rpt1').feature('sec2').feature('sec2').feature('geom1').set('children', {'sph1' 'off';  ...
'blk1' 'off';  ...
'blk2' 'off';  ...
'blk5' 'off';  ...
'imp1' 'off';  ...
'mov1' 'off';  ...
'blk3' 'off';  ...
'co1' 'off';  ...
'fin' 'off'});
model.result.report('rpt1').feature('sec2').feature('sec3').label('Materials');
model.result.report('rpt1').feature('sec2').feature('sec3').feature('sec1').label('Material 1');
model.result.report('rpt1').feature('sec2').feature('sec3').feature('sec1').set('source', 'firstchild');
model.result.report('rpt1').feature('sec2').feature('sec3').feature('sec1').feature('mat1').set('includesettings', false);
model.result.report('rpt1').feature('sec2').feature('sec3').feature('sec1').feature('mat1').set('children', {'def' 'off' 'off'});
model.result.report('rpt1').feature('sec2').feature('sec3').feature('sec2').label('Material 2');
model.result.report('rpt1').feature('sec2').feature('sec3').feature('sec2').set('source', 'firstchild');
model.result.report('rpt1').feature('sec2').feature('sec3').feature('sec2').feature('mat1').set('noderef', 'mat2');
model.result.report('rpt1').feature('sec2').feature('sec3').feature('sec2').feature('mat1').set('includesettings', false);
model.result.report('rpt1').feature('sec2').feature('sec3').feature('sec2').feature('mat1').set('children', {'def' 'off' 'off'});
model.result.report('rpt1').feature('sec2').feature('sec3').feature('sec3').label('Material 3');
model.result.report('rpt1').feature('sec2').feature('sec3').feature('sec3').set('source', 'firstchild');
model.result.report('rpt1').feature('sec2').feature('sec3').feature('sec3').feature('mat1').set('noderef', 'mat3');
model.result.report('rpt1').feature('sec2').feature('sec3').feature('sec3').feature('mat1').set('includesettings', false);
model.result.report('rpt1').feature('sec2').feature('sec3').feature('sec3').feature('mat1').set('children', {'def' 'off' 'off'});
model.result.report('rpt1').feature('sec2').feature('sec4').label('Electric Currents');
model.result.report('rpt1').feature('sec2').feature('sec4').set('source', 'firstchild');
model.result.report('rpt1').feature('sec2').feature('sec4').feature('phys1').set('includeusedproducts', false);
model.result.report('rpt1').feature('sec2').feature('sec4').feature('phys1').set('includeselection', false);
model.result.report('rpt1').feature('sec2').feature('sec4').feature('phys1').set('includesettings', false);
model.result.report('rpt1').feature('sec2').feature('sec4').feature('phys1').set('includefeaturetable', true);
model.result.report('rpt1').feature('sec2').feature('sec4').feature('phys1').set('children', {'cucn1' 'off' 'off' 'off';  ...
'ein1' 'off' 'off' 'off';  ...
'init1' 'off' 'off' 'off';  ...
'pot1' 'off' 'off' 'off';  ...
'dimp1' 'off' 'off' 'off';  ...
'gnd1' 'off' 'off' 'off'});
model.result.report('rpt1').feature('sec2').feature('sec4').feature('phys1').set('includevariables', false);
model.result.report('rpt1').feature('sec2').feature('sec4').feature('phys1').set('includeshapefunctions', false);
model.result.report('rpt1').feature('sec2').feature('sec4').feature('phys1').set('includeweakexpressions', false);
model.result.report('rpt1').feature('sec2').feature('sec4').feature('phys1').set('includeconstraints', false);
model.result.report('rpt1').feature('sec2').feature('sec5').label('Electric Currents 2');
model.result.report('rpt1').feature('sec2').feature('sec5').set('source', 'firstchild');
model.result.report('rpt1').feature('sec2').feature('sec5').feature('phys1').set('includeusedproducts', false);
model.result.report('rpt1').feature('sec2').feature('sec5').feature('phys1').set('noderef', 'ec2');
model.result.report('rpt1').feature('sec2').feature('sec5').feature('phys1').set('includeselection', false);
model.result.report('rpt1').feature('sec2').feature('sec5').feature('phys1').set('includesettings', false);
model.result.report('rpt1').feature('sec2').feature('sec5').feature('phys1').set('includefeaturetable', true);
model.result.report('rpt1').feature('sec2').feature('sec5').feature('phys1').set('children', {'cucn1' 'off' 'off' 'off';  ...
'ein1' 'off' 'off' 'off';  ...
'init1' 'off' 'off' 'off';  ...
'dimp1' 'off' 'off' 'off'});
model.result.report('rpt1').feature('sec2').feature('sec5').feature('phys1').set('includevariables', false);
model.result.report('rpt1').feature('sec2').feature('sec5').feature('phys1').set('includeshapefunctions', false);
model.result.report('rpt1').feature('sec2').feature('sec5').feature('phys1').set('includeweakexpressions', false);
model.result.report('rpt1').feature('sec2').feature('sec5').feature('phys1').set('includeconstraints', false);
model.result.report('rpt1').feature('sec2').feature('sec6').label('Weak Form Boundary PDE');
model.result.report('rpt1').feature('sec2').feature('sec6').set('source', 'firstchild');
model.result.report('rpt1').feature('sec2').feature('sec6').feature('phys1').set('includeusedproducts', false);
model.result.report('rpt1').feature('sec2').feature('sec6').feature('phys1').set('noderef', 'wb');
model.result.report('rpt1').feature('sec2').feature('sec6').feature('phys1').set('includeselection', false);
model.result.report('rpt1').feature('sec2').feature('sec6').feature('phys1').set('includesettings', false);
model.result.report('rpt1').feature('sec2').feature('sec6').feature('phys1').set('includefeaturetable', true);
model.result.report('rpt1').feature('sec2').feature('sec6').feature('phys1').set('children', {'wfeq1' 'off' 'off' 'off'; 'init1' 'off' 'off' 'off'});
model.result.report('rpt1').feature('sec2').feature('sec6').feature('phys1').set('includevariables', false);
model.result.report('rpt1').feature('sec2').feature('sec6').feature('phys1').set('includeshapefunctions', false);
model.result.report('rpt1').feature('sec2').feature('sec6').feature('phys1').set('includeweakexpressions', false);
model.result.report('rpt1').feature('sec2').feature('sec6').feature('phys1').set('includeconstraints', false);
model.result.report('rpt1').feature('sec2').feature('sec7').label('Mesh 1');
model.result.report('rpt1').feature('sec2').feature('sec7').set('source', 'firstchild');
model.result.report('rpt1').feature('sec2').feature('sec7').feature('mesh1').set('children', {'size' 'off'; 'ftet2' 'off'; 'ftet2' 'off'; 'ftet1' 'off'; 'ftet3' 'off'});
model.result.report('rpt1').feature('sec3').label('Study 1');
model.result.report('rpt1').feature('sec4').label('Results');
model.result.report('rpt1').feature('sec4').feature('sec1').label('Data Sets');
model.result.report('rpt1').feature('sec4').feature('sec1').feature('sec1').label('Study 2 MeshSize/Solution 6');
model.result.report('rpt1').feature('sec4').feature('sec1').feature('sec1').set('source', 'firstchild');
model.result.report('rpt1').feature('sec4').feature('sec1').feature('sec1').feature('dset1').set('includesettings', false);
model.result.report('rpt1').feature('sec4').feature('sec1').feature('sec2').label('Study 2 MeshSize/Solution 6.1');
model.result.report('rpt1').feature('sec4').feature('sec1').feature('sec2').set('source', 'firstchild');
model.result.report('rpt1').feature('sec4').feature('sec1').feature('sec2').feature('dset1').set('includesettings', false);
model.result.report('rpt1').feature('sec4').feature('sec2').label('Tables');
model.result.report('rpt1').feature('sec4').feature('sec2').feature('sec1').label('Evaluation 3D');
model.result.report('rpt1').feature('sec4').feature('sec2').feature('sec1').set('source', 'firstchild');
model.result.report('rpt1').feature('sec4').feature('sec3').label('Plot Groups');
model.result.report('rpt1').feature('sec4').feature('sec3').feature('sec1').label('Influence-MeshSize');
model.result.report('rpt1').feature('sec4').feature('sec3').feature('sec1').set('source', 'firstchild');
model.result.report('rpt1').feature('sec4').feature('sec3').feature('sec1').feature('pg1').set('noderef', 'pg14');
model.result.report('rpt1').feature('sec4').feature('sec3').feature('sec2').label('Influence-MeshSize 1');
model.result.report('rpt1').feature('sec4').feature('sec3').feature('sec2').set('source', 'firstchild');
model.result.report('rpt1').feature('sec4').feature('sec3').feature('sec2').feature('pg1').set('noderef', 'pg14');
model.result.report('rpt1').feature('sec4').feature('sec3').feature('sec3').label('Influence-MeshSize 1');
model.result.report('rpt1').feature('sec4').feature('sec3').feature('sec3').set('source', 'firstchild');
model.result.report('rpt1').feature('sec4').feature('sec3').feature('sec3').feature('pg1').set('noderef', 'pg14');
model.result.report('rpt1').feature('sec4').feature('sec3').feature('sec4').label('Influence-MeshSize 1');
model.result.report('rpt1').feature('sec4').feature('sec3').feature('sec4').set('source', 'firstchild');
model.result.report('rpt1').feature('sec4').feature('sec3').feature('sec4').feature('pg1').set('noderef', 'pg14');
model.result.report('rpt1').feature('sec4').feature('sec3').feature('sec5').label('Influence-MeshSize 1');
model.result.report('rpt1').feature('sec4').feature('sec3').feature('sec5').set('source', 'firstchild');
model.result.report('rpt1').feature('sec4').feature('sec3').feature('sec5').feature('pg1').set('noderef', 'pg14');

model.label('Microtrap_diff_2_output_solved_LR.mph');

model.study.remove('std1');
model.study.remove('std2');
model.study.create('std4');
model.study('std4').create('time', 'Transient');
model.study('std4').feature('time').activate('ec', true);
model.study('std4').feature('time').activate('ec2', true);
model.study('std4').feature('time').activate('wb', true);
model.study('std4').feature('time').set('tlist', '10^{range(log10(1.0e-8),(log10(tpulse[1/s])+8)/100,log10(tpulse[1/s]))}');
model.study('std4').feature('time').set('usertol', true);
model.study('std4').feature('time').set('rtol', '0.0001');
model.study('std4').create('param', 'Parametric');
model.study('std4').feature('param').setIndex('pname', 'sigma_e', 0);
model.study('std4').feature('param').setIndex('plistarr', '', 0);
model.study('std4').feature('param').setIndex('punit', 'S/m', 0);
model.study('std4').feature('param').setIndex('pname', 'sigma_e', 0);
model.study('std4').feature('param').setIndex('plistarr', '', 0);
model.study('std4').feature('param').setIndex('punit', 'S/m', 0);
model.study('std4').feature('param').setIndex('pname', 'Uapp', 0);
model.study('std4').feature('param').setIndex('plistarr', '3.595, 7.3644, 11.047', 0);

model.sol.create('sol23');
model.sol('sol23').study('std4');

model.study('std4').feature('time').set('notlistsolnum', 1);
model.study('std4').feature('time').set('notsolnum', '1');
model.study('std4').feature('time').set('listsolnum', 1);
model.study('std4').feature('time').set('solnum', '1');

model.sol('sol23').create('st1', 'StudyStep');
model.sol('sol23').feature('st1').set('study', 'std4');
model.sol('sol23').feature('st1').set('studystep', 'time');
model.sol('sol23').create('v1', 'Variables');
model.sol('sol23').feature('v1').set('control', 'time');
model.sol('sol23').create('t1', 'Time');
model.sol('sol23').feature('t1').set('tlist', '10^{range(log10(1.0e-8),(log10(tpulse[1/s])+8)/100,log10(tpulse[1/s]))}');
model.sol('sol23').feature('t1').set('plot', 'off');
model.sol('sol23').feature('t1').set('plotgroup', 'pg18');
model.sol('sol23').feature('t1').set('plotfreq', 'tout');
model.sol('sol23').feature('t1').set('probesel', 'all');
model.sol('sol23').feature('t1').set('probes', {});
model.sol('sol23').feature('t1').set('probefreq', 'tsteps');
model.sol('sol23').feature('t1').set('atolglobalvaluemethod', 'factor');
model.sol('sol23').feature('t1').set('endtimeinterpolation', true);
model.sol('sol23').feature('t1').set('control', 'time');
model.sol('sol23').feature('t1').create('se1', 'Segregated');
model.sol('sol23').feature('t1').feature('se1').feature.remove('ssDef');
model.sol('sol23').feature('t1').feature('se1').create('ss1', 'SegregatedStep');
model.sol('sol23').feature('t1').feature('se1').feature('ss1').set('segvar', {'comp1_V'});
model.sol('sol23').feature('t1').create('i1', 'Iterative');
model.sol('sol23').feature('t1').feature('i1').set('linsolver', 'cg');
model.sol('sol23').feature('t1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol23').feature('t1').feature('i1').feature('mg1').set('prefun', 'amg');
model.sol('sol23').feature('t1').feature('se1').feature('ss1').set('linsolver', 'i1');
model.sol('sol23').feature('t1').feature('se1').feature('ss1').label('Electric Currents');
model.sol('sol23').feature('t1').feature('se1').create('ss2', 'SegregatedStep');
model.sol('sol23').feature('t1').feature('se1').feature('ss2').set('segvar', {'comp1_V2'});
model.sol('sol23').feature('t1').create('i2', 'Iterative');
model.sol('sol23').feature('t1').feature('i2').set('linsolver', 'cg');
model.sol('sol23').feature('t1').feature('i2').create('mg1', 'Multigrid');
model.sol('sol23').feature('t1').feature('i2').feature('mg1').set('prefun', 'amg');
model.sol('sol23').feature('t1').feature('se1').feature('ss2').set('linsolver', 'i2');
model.sol('sol23').feature('t1').feature('se1').feature('ss2').label('Electric Currents 2');
model.sol('sol23').feature('t1').feature('se1').create('ss3', 'SegregatedStep');
model.sol('sol23').feature('t1').feature('se1').feature('ss3').set('segvar', {'comp1_u'});
model.sol('sol23').feature('t1').feature('se1').feature('ss3').set('linsolver', 'dDef');
model.sol('sol23').feature('t1').feature('se1').feature('ss3').label('Weak Form Boundary PDE');
model.sol('sol23').feature('t1').feature.remove('fcDef');
model.sol('sol23').attach('std4');

model.batch.create('p5', 'Parametric');
model.batch('p5').study('std4');
model.batch('p5').create('so1', 'Solutionseq');
model.batch('p5').feature('so1').set('seq', 'sol23');
model.batch('p5').feature('so1').set('store', 'on');
model.batch('p5').feature('so1').set('clear', 'on');
model.batch('p5').feature('so1').set('psol', 'none');
model.batch('p5').set('pname', {'Uapp'});
model.batch('p5').set('plistarr', {'3.595, 7.3644, 11.047'});
model.batch('p5').set('sweeptype', 'sparse');
model.batch('p5').set('probesel', 'all');
model.batch('p5').set('probes', {});
model.batch('p5').set('plot', 'off');
model.batch('p5').set('err', 'on');
model.batch('p5').attach('std4');
model.batch('p5').set('control', 'param');

model.sol('sol23').feature('t1').set('tstepsbdf', 'strict');
model.sol('sol23').feature('t1').set('initialstepbdfactive', true);
model.sol('sol23').feature('t1').set('initialstepbdf', '1e-9');

model.study('std4').feature('time').set('timeadaption', true);

model.component('comp1').material('mat4').selection.set([]);
model.component('comp1').material('mat1').selection.set([1 3]);

model.study('std4').setGenPlots(false);
model.study('std4').setGenConv(false);

model.sol.create('sol24');
model.sol('sol24').label('Refined Mesh Solution 2');
model.sol('sol24').study('std4');
model.sol('sol24').setClusterStorage('all');
model.sol('sol23').feature('t1').feature('taDef').set('tadapsol', 'sol24');
model.sol.create('sol25');
model.sol('sol25').study('std4');
model.sol('sol25').label('Parametric Solutions 1');

model.batch('p5').feature('so1').set('psol', 'sol25');

model.sol('sol23').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol23').feature('t1').feature('dDef').active(true);
model.sol('sol23').feature('t1').feature.remove('se1');
model.sol('sol23').feature('t1').feature.remove('i1');
model.sol('sol23').feature('t1').feature.remove('i2');

model.batch('p5').run;

model.study('std4').label('Study 4 Trap=ExtMedium');

model.result('pg18').label('Pore density');
model.result.remove('pg20');
model.result('pg18').feature('lngr2').label('With trap');
model.result('pg18').feature.duplicate('lngr3', 'lngr2');
model.result('pg18').feature('lngr3').label('Without trap');
model.result('pg18').feature('lngr3').set('data', 'dset12');
model.result('pg18').feature('lngr2').set('linecolor', 'cycle');
model.result('pg18').feature('lngr2').set('linewidth', 4);
model.result('pg18').feature('lngr2').set('linestyle', 'solid');
model.result('pg18').feature('lngr3').set('linecolor', 'cyclereset');
model.result('pg18').run;
model.result('pg18').feature('lngr3').set('linewidth', 3);

model.label('Microtrap_diff_2_output_solved_LR_trapVSnotrap.mph');

model.param.set('rp', '1e-9[m]');

model.component('comp1').variable('var1').set('ITV', 'V-V2');
model.component('comp1').variable.duplicate('var4', 'var2');
model.component('comp1').variable('var2').active(false);
model.component('comp1').variable('var4').remove('K0');
model.component('comp1').variable('var4').remove('K1');
model.component('comp1').variable('var4').remove('K2');
model.component('comp1').variable('var4').remove('inter');
model.component('comp1').variable('var4').remove('sigma_ep');
model.component('comp1').variable('var4').remove('vm');
model.component('comp1').variable('var4').rename('sigma_ep_new', 'sigma_ep');

model.batch('p5').run;

model.component('comp1').material('mat4').selection.set([3]);
model.component('comp1').material('mat1').selection.set([1]);

model.batch('p4').run;

model.result('pg18').run;

model.component('comp1').material('mat4').selection.set([]);
model.component('comp1').material('mat1').selection.set([1 3]);

model.batch('p5').run;

model.result('pg18').run;
model.result('pg18').feature('lngr3').setIndex('looplevelinput', 'last', 0);
model.result('pg18').run;

model.label('Microtrap_diff_2_output_solved_LR_trapVSnotrap_updateEPmodel.mph');

out = model;
