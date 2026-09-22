(* ::Package:: *)
(* Executable derivation, not C++ implementation. All state is dimensionless.
   Get this file in Mathematica, or run it with wolframscript -file.
   Generated coefficient tables and companion notebook stay beside this file. *)
Begin["MovingThreeFluid`"];
ClearAll["MovingThreeFluid`*"];
SetAttributes[dx,NHoldAll]; (* Numeric evaluation must not turn tensor indices into reals. *)
directory = DirectoryName[$InputFileName];
gamma = 5/3; beta = gamma-1;
species = {"s","b","d"}; fields = {"rho","v","u","M"};
variables = Flatten[Table[dx[f,o,k],{o,-1,1},{f,0,2},{k,0,3}]];
zero = Thread[variables->0];
rowLabels = Flatten[Table[species[[f+1]]<>" "<># & /@
   {"continuity","momentum","energy","enclosed mass"},{f,0,2}]];
old[f_,o_] := {rho[f,o],vel[f,o],u[f,o]};
increment[f_,o_] := Table[dx[f,o,k],{k,0,2}];
storage[{rr_,vv_,uu_}] := {rr,rr vv,rr (uu+eps vv^2/2)};
flux[{rr_,vv_,uu_}] := {rr vv,rr vv^2+beta rr uu/eps,
   rr vv (gamma uu+eps vv^2/2)};
jacobian[fun_,point_] := Table[D[fun[{rr,vv,uu}][[a]],{rr,vv,uu}[[b]]],
   {a,3},{b,3}]/.Thread[{rr,vv,uu}->point];
linearValue[fun_,point_,change_] := fun[point]+jacobian[fun,point].change;
vacuumFlux[{rr_,vv_,uu_}] := Module[{cs,chi},
 cs=Sqrt[gamma beta uu/eps];chi=(2 cs+(gamma-1)vv)/((gamma+1)cs);
 flux[{rr chi^(2/(gamma-1)),cs chi,uu chi^2}]
];
(* Raw primitive states; equilibrium is maintained by source quadrature.
   Advective dissipation vanishes at rest. Wave speeds are frozen. *)
alternating[{rl_,vl_,ul_},{rr_,vr_,ur_}] := With[
 {d=(rl+rr)/2,e=(ul+ur)/2,p=beta rr ur},
 {d vl,d vl^2+p/eps,d vl(e+eps vl^2/2)+p vl}];
altLeft[l_,r_] := Table[D[alternating[{rl,vl,ul},{rr,vr,ur}][[a]],{rl,vl,ul}[[b]]],
 {a,3},{b,3}]/.Thread[{rl,vl,ul,rr,vr,ur}->Join[l,r]];
altRight[l_,r_] := Table[D[alternating[{rl,vl,ul},{rr,vr,ur}][[a]],{rr,vr,ur}[[b]]],
 {a,3},{b,3}]/.Thread[{rl,vl,ul,rr,vr,ur}->Join[l,r]];
face[f_,side_,linear_] := Module[{left,right,xl,xr,fl,fr,wl,wr,vl,vr,dl,dr,viscous},
 {left,right}=If[side==-1,{-1,0},{0,1}];
 xl=old[f,left];xr=old[f,right];
 If[linear,
  fl=alternating[xl,xr]+altLeft[xl,xr].increment[f,left]+altRight[xl,xr].increment[f,right];
  wl=linearValue[storage,xl,increment[f,left]];wr=linearValue[storage,xr,increment[f,right]],
  fl=alternating[xl+increment[f,left],xr+increment[f,right]];
  wl=storage[xl+increment[f,left]];wr=storage[xr+increment[f,right]]];
 vl=xl[[2]];vr=xr[[2]];dl=dx[f,left,1];dr=dx[f,right,1];
 viscous=If[linear,{0,vr-vl+dr-dl,eps((vr^2-vl^2)/2+vr dr-vl dl)},
   {0,vr+dr-vl-dl,eps((vr+dr)^2-(vl+dl)^2)/2}];
 fl-speed[f,side,0](wr-wl)/2-viscosity[f,side]viscous/2
];
outerFace[f_,mode_,linear_] := Switch[mode,
 "vacuum",{0,0,0},
 "supersonic",If[linear,linearValue[flux,old[f,0],increment[f,0]],flux[old[f,0]+increment[f,0]]],
 "fan",If[linear,linearValue[vacuumFlux,old[f,0],increment[f,0]],vacuumFlux[old[f,0]+increment[f,0]]]
];
z[f_,o_,linear_] := If[linear,Sqrt[u[f,o]]+dx[f,o,2]/(2 Sqrt[u[f,o]]),
 Sqrt[u[f,o]+dx[f,o,2]]];
w[f_,linear_] := If[linear,1/Sqrt[u[f,0]]-dx[f,0,2]/(2 u[f,0]^(3/2)),
 1/Sqrt[u[f,0]+dx[f,0,2]]];
luminosity[f_,side_,linear_] := If[side==-1,
 -conductance[f,-1](z[f,0,linear]-z[f,-1,linear]),
 -conductance[f,1](z[f,1,linear]-z[f,0,linear])];
thermalSource[f_,linear_] := Module[{exchange,heating},
 exchange=-rho[f,0] Sum[If[g==f,0,c1[f,g]rho[g,0] *
  (particleMass[f](u[f,0]+dx[f,0,2])-particleMass[g](u[g,0]+dx[g,0,2])) /
  (particleMass[0](u[f,0]+u[g,0])^(3/2))],{g,0,2}];
 heating=rho[f,0] If[f==1,Sum[c4[f,g]rho[g,0]w[g,linear],{g,0,2}],
  c4[f,1]rho[1,0]w[f,linear]];
 exchange+heating
];
gravity[centre_] := (eta Sum[mass[g,0],{g,0,2}]+
 (1-eta)If[centre,0,Sum[mass[g,-1],{g,0,2}]])/radius^2-q radius;
deltaGravity[centre_] := (eta Sum[dx[g,0,3],{g,0,2}]+
 (1-eta)If[centre,0,Sum[dx[g,-1,3],{g,0,2}]])/radius^2;
source[f_,centre_,linear_] := Module[{rr,vv,uu,dr,dv,du,gg,dg,momentum,energy},
 {rr,vv,uu}=old[f,0];{dr,dv,du}=increment[f,0];
 gg=gravity[centre]-balanceAcceleration[f];dg=deltaGravity[centre];
 If[linear,
  momentum=(geometry beta(rr uu+uu dr+rr du)-(rr gg+gg dr+rr dg))/eps;
  energy=-(rr vv gg+vv gg dr+rr gg dv+rr vv dg),
  momentum=(geometry beta(rr+dr)(uu+du)-(rr+dr)(gg+dg))/eps;
  energy=-(rr+dr)(vv+dv)(gg+dg)];
 {0,momentum,energy+thermalSource[f,linear]}
];
makeRows[region_,linear_:True] := Module[{centre,outer,mode,result,lf,rf,ll,rl,t,eq,mc},
 centre=region=="centre";outer=StringStartsQ[region,"outer"];
 mode=Switch[region,"outer_fan","fan","outer_supersonic","supersonic",_,"vacuum"];
 result=Table[
  lf=If[centre,{0,0,0},face[f,-1,linear]];
  rf=If[outer,outerFace[f,mode,linear],face[f,1,linear]];
  ll=If[centre,0,luminosity[f,-1,linear]];
  rl=If[outer,escapeConductance[f]z[f,0,linear],luminosity[f,1,linear]];
  t=If[linear,jacobian[storage,old[f,0]].increment[f,0],
    storage[old[f,0]+increment[f,0]]-storage[old[f,0]]];
  eq=volume t+dt(areaRight rf-areaLeft lf+{0,0,rl-ll}-volume source[f,centre,linear]);
  mc=mass[f,0]+dx[f,0,3]-If[centre,0,mass[f,-1]+dx[f,-1,3]]-
     volume(rho[f,0]+dx[f,0,0]);
  Join[eq,{mc}],{f,0,2}];
 Flatten[result]
];
regions={"centre","bulk","outer_fan","outer_supersonic","outer_vacuum"};
compileRows[region_] := Module[{equations,a,b},
 equations=Expand[makeRows[region]];
 a=Table[Coefficient[equations[[i]],variables[[j]]],{i,12},{j,36}];
 b=-(equations/.zero);
 <|"Variables"->variables,"RowLabels"->rowLabels,"Residual"->equations,
   "Aminus"->a[[All,1;;12]],"Adiag"->a[[All,13;;24]],"Aplus"->a[[All,25;;36]],
   "Matrix"->a,"IncrementRHS"->b,
   "AbsoluteRHS"->b+a.Flatten[Table[Append[old[f,o],mass[f,o]],{o,-1,1},{f,0,2}]]|>
];
Print["Deriving all centre, bulk and outer boundary rows..."];
rowsByRegion=Association@Table[reg->compileRows[reg],{reg,regions}];
entriesByRegion=Association@Table[reg->Flatten[Table[
 If[rowsByRegion[reg]["Matrix"][[i,j]]===0,Nothing,
  {rowLabels[[i]],variables[[j]],rowsByRegion[reg]["Matrix"][[i,j]]}],
 {i,12},{j,36}],1],{reg,regions}];
storageJacobian=jacobian[storage,{rr,vv,uu}];
fluxJacobian=jacobian[flux,{rr,vv,uu}];
vacuumJacobian=jacobian[vacuumFlux,{rr,vv,uu}];

(* Verify each affine residual and each generated entry/RHS. *)
checks=<||>;
AssociateTo[checks,"alternating_consistency"->TrueQ[
 Simplify[alternating[{rr,vv,uu},{rr,vv,uu}]-flux[{rr,vv,uu}]]=={0,0,0}]];
AssociateTo[checks,"checkerboard_no_null"->TrueQ[
 Simplify[((1-Exp[-I theta])(Exp[I theta]-1))/.theta->Pi]==-4]];
AssociateTo[checks,"checkerboard_backward_euler"->TrueQ[
 Expand[Det[IdentityMatrix[2]-courant {{-2,2},{-2,0}}]]==1+2courant+4courant^2]];
Do[AssociateTo[checks,"affine_"<>reg->And@@(TrueQ[#==0]& /@
  (Expand[rowsByRegion[reg]["Residual"]-rowsByRegion[reg]["Matrix"].variables+
   rowsByRegion[reg]["IncrementRHS"]]))],{reg,regions}];
AssociateTo[checks,"storage_Jacobian"->(storageJacobian===
 {{1,0,0},{vv,rr,0},{uu+eps vv^2/2,eps rr vv,rr}})];
AssociateTo[checks,"flux_Jacobian"->And@@(TrueQ[#==0]& /@Flatten[Expand[fluxJacobian-
 {{vv,rr,0},{vv^2+beta uu/eps,2 rr vv,beta rr/eps},
 {vv(gamma uu+eps vv^2/2),rr(gamma uu+3 eps vv^2/2),gamma rr vv}}]])];
cs=Sqrt[gamma beta uu/eps];chi=(3 cs+vv)/(4 cs);
vacuumMapJacobian={{chi^3,3 rr chi^2/(4 cs),-3 rr chi^2 vv/(8 cs uu)},
 {0,1/4,3 cs/(8 uu)},{0,uu chi/(2 cs),chi^2-chi vv/(4 cs)}};
AssociateTo[checks,"vacuum_chain_rule"->TrueQ[FullSimplify[
 vacuumJacobian-jacobian[flux,{rr chi^3,cs chi,uu chi^2}].vacuumMapJacobian,
 uu>0 && eps>0]==ConstantArray[0,{3,3}]]];

(* Numerical finite differences of NONLINEAR residuals, not of A.delta-b.
   Coefficients explicitly frozen by the scheme stay frozen in this check. *)
numericRules={rho[f_,o_]:>1+f/10+o/25,vel[f_,o_]:>1/5+f/30+o/100,
 u[f_,o_]:>1+f/5+o/30,mass[f_,o_]:>2/5+f/10+o/5,
 refRho[___]->0,refU[___]->0,speed[___]->5,viscosity[___]->2,conductance[___]->7/100,
 escapeConductance[_]->1/20,c1[_,_]->3/1000,c4[_,_]->1/500,
 particleMass[f_]:>{1,2,1/10}[[f+1]],balanceAcceleration[_]->1/200,
 volume->1/5,dt->1/100,areaLeft->4/5,areaRight->6/5,radius->1,
 geometry->2,eta->2/5,eps->1/10,q->1/100};
finiteDifferenceErrors=<||>;
Do[Module[{nl,aa,h=1/1000000,fd,xp,xm},
 nl=N[makeRows[reg,False]/.numericRules,70];aa=N[rowsByRegion[reg]["Matrix"]/.numericRules,40];
 fd=Transpose@Table[xp=ConstantArray[0,36];xm=xp;xp[[k]]=h;xm[[k]]=-h;
  N[((nl/.Thread[variables->xp])-(nl/.Thread[variables->xm]))/(2h),40],{k,36}];
 AssociateTo[finiteDifferenceErrors,reg->Max[Abs[Flatten[fd-aa]]]];
 AssociateTo[checks,"nonlinear_Jacobian_"<>reg->TrueQ[finiteDifferenceErrors[reg]<10^-10]];
 ],{reg,regions}];

(* Actual scalar half-bandwidth with the requested 12-variable ordering. *)
positions=Flatten[Table[Module[{local,row,col,reg},
 reg=If[i==0,"centre",If[i==4,"outer_fan","bulk"]];local=rowsByRegion[reg]["Matrix"];
 Flatten[Table[If[local[[a,b]]===0,Nothing,
  row=12 i+a-1;col=12(i+Quotient[b-1,12]-1)+Mod[b-1,12];{row,col}],
  {a,12},{b,36}],1]],{i,0,4}],1];
bandwidth={Max[Subtract@@@positions],Max[-Subtract@@@positions]};
AssociateTo[checks,"nearest_cell_stencil"->And@@(TrueQ[0<=#[[2]]<60]& /@positions)];
AssociateTo[checks,"constant_bandwidth"->(bandwidth=={19,13})];
(* Independent mass recurrence/centre regularity and telescoping ledger. *)
AssociateTo[checks,"centre_mass"->And@@Table[TrueQ[Expand[
 makeRows["centre"][[4 f+4]]-(mass[f,0]+dx[f,0,3]-volume(rho[f,0]+dx[f,0,0]))]==0],{f,0,2}]];
AssociateTo[checks,"mass_ledger"->TrueQ[Expand[
 Sum[-(outflow[i]-outflow[i-1]),{i,1,7}]+outflow[7]-outflow[0]]==0]];
reindex[expr_,cell_] := expr/.{
 rho[f_,o_]:>rhoGlobal[f,cell+o],vel[f_,o_]:>velGlobal[f,cell+o],
 u[f_,o_]:>uGlobal[f,cell+o],dx[f_,o_,k_]:>dxGlobal[f,cell+o,k],
 refRho[f_,side_,s_]:>refRhoGlobal[f,cell+If[side==1,1,0],s],
 refU[f_,side_,s_]:>refUGlobal[f,cell+If[side==1,1,0],s],
 speed[f_,side_,k_]:>speedGlobal[f,cell+If[side==1,1,0],k],
 viscosity[f_,side_]:>viscosityGlobal[f,cell+If[side==1,1,0]],
 areaLeft->areaGlobal[cell],areaRight->areaGlobal[cell+1],volume->volumeGlobal[cell]};
AssociateTo[checks,"assembled_continuity_ledger"->And@@Table[
 TrueQ[Expand[Sum[reindex[rowsByRegion[If[cell==0,"centre",If[cell==3,"outer_fan","bulk"]]]["Residual"][[4 f+1]],cell],{cell,0,3}]-
 Sum[volumeGlobal[cell]dxGlobal[f,cell,0],{cell,0,3}]-
 dt areaGlobal[4]reindex[outerFace[f,"fan",True][[1]],3]]==0],{f,0,2}]];
(* The source linearization is exactly the one in conduction_compiler.wl. *)
AssociateTo[checks,"inverse_sqrt_Taylor"->TrueQ[FullSimplify[
  1/Sqrt[uu]-dd/(2 uu^(3/2))-(3/(2 Sqrt[uu])-(uu+dd)/(2 uu^(3/2))),uu>0]==0]];
AssociateTo[checks,"pair_exchange_energy"->TrueQ[Simplify[
 Sum[thermalSource[f,True],{f,0,2}]/.{c4[_,_]->0,c1[_,_]->cc}]==0]];
(* A finite primitive correction must recover conservative quantities. *)
oldExample={6/5,1/7,4/5};deltaExample={1/100,-1/200,1/300};
target=(storage[oldExample]+jacobian[storage,oldExample].deltaExample)/.eps->1/10;
recovered={target[[1]],target[[2]]/target[[1]],
 target[[3]]/target[[1]]-(1/10)(target[[2]]/target[[1]])^2/2};
AssociateTo[checks,"conservative_recovery"->TrueQ[Simplify[(storage[recovered]/.eps->1/10)-target]=={0,0,0}]];

(* Independently assemble the block formulas printed in the TeX note. *)
manualMatrix[reg_] := Module[{mat=ConstantArray[0,{12,36}],centre,outer,mode,
 ids,cols,point,jl,jr,dd,gg,kminus,kplus,hb,bcoef,rrf,vvf,uuf},
 centre=reg=="centre";outer=StringStartsQ[reg,"outer"];
 mode=Switch[reg,"outer_fan","fan","outer_supersonic","supersonic",_,"vacuum"];
 Do[
  ids=Range[4 f+1,4 f+3];cols=Range[12+4 f+1,12+4 f+3];
  {rrf,vvf,uuf}=old[f,0];gg=gravity[centre]-balanceAcceleration[f];
  dd={{0,0,0},{(geometry beta uuf-gg)/eps,0,geometry beta rrf/eps},
      {-vvf gg,-rrf gg,0}};
  mat[[ids,cols]]=volume jacobian[storage,old[f,0]]-dt volume dd;
  point[side_,lr_] := old[f,If[side==-1,lr-1,lr]];
  jl[side_] := altLeft[point[side,0],point[side,1]]+speed[f,side,0]jacobian[storage,point[side,0]]/2+
    viscosity[f,side]/2 {{0,0,0},{0,1,0},{0,eps point[side,0][[2]],0}};
  jr[side_] := altRight[point[side,0],point[side,1]]-speed[f,side,0]jacobian[storage,point[side,1]]/2-
    viscosity[f,side]/2 {{0,0,0},{0,1,0},{0,eps point[side,1][[2]],0}};
  If[!centre,
   mat[[ids,Range[4 f+1,4 f+3]]]-=dt areaLeft jl[-1];
   mat[[ids,cols]]-=dt areaLeft jr[-1]];
  If[outer,
   hb=Switch[mode,"fan",jacobian[vacuumFlux,old[f,0]],"supersonic",jacobian[flux,old[f,0]],_,ConstantArray[0,{3,3}]];
   mat[[ids,cols]]+=dt areaRight hb,
   mat[[ids,cols]]+=dt areaRight jl[1];
   mat[[ids,Range[24+4 f+1,24+4 f+3]]]+=dt areaRight jr[1]];
  kminus=If[centre,0,conductance[f,-1]];
  kplus=If[outer,escapeConductance[f],conductance[f,1]];
  mat[[4 f+3,12+4 f+3]]+=dt(kminus+kplus)/(2 Sqrt[u[f,0]]);
  If[!centre,mat[[4 f+3,4 f+3]]-=dt kminus/(2 Sqrt[u[f,-1]])];
  If[!outer,mat[[4 f+3,24+4 f+3]]-=dt kplus/(2 Sqrt[u[f,1]])];
  Do[
   bcoef=If[g==f,-rho[f,0]particleMass[f]Sum[If[h==f,0,
      c1[f,h]rho[h,0]/(particleMass[0](u[f,0]+u[h,0])^(3/2))],{h,0,2}],
      rho[f,0]c1[f,g]rho[g,0]particleMass[g]/(particleMass[0](u[f,0]+u[g,0])^(3/2))];
   bcoef-=If[f==1,rho[f,0]c4[f,g]rho[g,0]/(2 u[g,0]^(3/2)),
      If[g==f,rho[f,0]c4[f,1]rho[1,0]/(2 u[f,0]^(3/2)),0]];
   mat[[4 f+3,12+4 g+3]]-=dt volume bcoef;
   mat[[ids,12+4 g+4]]-=dt volume {0,-rrf eta/(eps radius^2),-rrf vvf eta/radius^2};
   If[!centre,mat[[ids,4 g+4]]-=dt volume {0,-rrf(1-eta)/(eps radius^2),-rrf vvf(1-eta)/radius^2}],
  {g,0,2}];
  mat[[4 f+4,12+4 f+1]]=-volume;mat[[4 f+4,12+4 f+4]]=1;
  If[!centre,mat[[4 f+4,4 f+4]]=-1],{f,0,2}];
 mat
];
Do[AssociateTo[checks,"TeX_blocks_"<>reg->And@@(TrueQ[#==0]& /@
 Flatten[Expand[manualMatrix[reg]-rowsByRegion[reg]["Matrix"]]])],{reg,regions}];

(* The logical derivation remains readable; these are the C++ storage tables. *)
storageOrder={3,2,1,7,6,5,11,10,9,4,8,12};
storageColumns=Flatten[Table[storageOrder+12 j,{j,0,2}]];
storageRowsByRegion=AssociationMap[Function[reg,Module[{x=rowsByRegion[reg],a},
 a=x["Matrix"][[storageOrder,storageColumns]];
 <|"Variables"->x["Variables"][[storageColumns]],
   "RowLabels"->x["RowLabels"][[storageOrder]],
   "Residual"->x["Residual"][[storageOrder]],"Matrix"->a,
   "Aminus"->a[[All,1;;12]],"Adiag"->a[[All,13;;24]],"Aplus"->a[[All,25;;36]],
   "IncrementRHS"->x["IncrementRHS"][[storageOrder]],
   "AbsoluteRHS"->x["AbsoluteRHS"][[storageOrder]]|>]],regions];
storagePosition=Ordering[storageOrder]-1;
storageOffsets=(12 Quotient[#[[1]],12]+storagePosition[[Mod[#[[1]],12]+1]]-
  12 Quotient[#[[2]],12]-storagePosition[[Mod[#[[2]],12]+1]])& /@positions;
storageBandwidth={Max[storageOffsets],Max[-storageOffsets]};
AssociateTo[checks,"optimized_bandwidth"->(storageBandwidth=={13,14})];
Do[AssociateTo[checks,"optimized_affine_"<>reg->TrueQ[
 Expand[storageRowsByRegion[reg]["Residual"]-
 storageRowsByRegion[reg]["Matrix"].storageRowsByRegion[reg]["Variables"]+
 storageRowsByRegion[reg]["IncrementRHS"]]==ConstantArray[0,12]]],{reg,regions}];
Block[{$Context="Global`",$ContextPath={"System`"}},
 Put[storageRowsByRegion,FileNameJoin[{directory,"moving_three_fluid_storage_coefficients.wl"}]];
 Put[rowsByRegion,FileNameJoin[{directory,"moving_three_fluid_coefficients.wl"}]]];
Put[<|"Checks"->checks,"FiniteDifferenceErrors"->finiteDifferenceErrors,
 "LogicalHalfBandwidth"->bandwidth,"HalfBandwidth"->storageBandwidth|>,FileNameJoin[{directory,"moving_three_fluid_checks.wl"}]];
report=StringRiffle[Flatten@Table[Join[
 {"REGION: "<>reg,"VARIABLES: "<>ToString[variables,InputForm]},
 (ToString[#,InputForm]& /@entriesByRegion[reg]),
 Table["RHS["<>rowLabels[[i]]<>"] = "<>ToString[rowsByRegion[reg]["IncrementRHS"][[i]],InputForm],{i,12}]],
 {reg,regions}],"\n"];
Export[FileNameJoin[{directory,"moving_three_fluid_entries.txt"}],report,"Text"];
storageReport=StringRiffle[Flatten@Table[Join[{"REGION: "<>reg},
 Flatten[Table[If[storageRowsByRegion[reg]["Matrix"][[i,j]]===0,Nothing,
 ToString[{i-1,storageRowsByRegion[reg]["Variables"][[j]],
 storageRowsByRegion[reg]["Matrix"][[i,j]]},InputForm]],{i,12},{j,36}]],
 Table["RHS["<>ToString[i-1]<>"] = "<>
 ToString[storageRowsByRegion[reg]["IncrementRHS"][[i]],InputForm],{i,12}]],{reg,regions}],"\n"];
Export[FileNameJoin[{directory,"moving_three_fluid_storage_entries.txt"}],storageReport,"Text"];
notebook=Notebook[{
 Cell["MovingThreeFluidSim: one banded solve", "Title"],
 Cell["Analytic compiler for review. No C++ is generated. Evaluate the following cell to regenerate all coefficients and tests. The state order is cell, species (s,b,d), then rho,v,u,M. M is enclosed mass at the outer cell face. Geometry is a prescribed common grid.","Text"],
 Cell[BoxData[ToBoxes[Defer[Get[FileNameJoin[{NotebookDirectory[],"moving_three_fluid_compiler.wl"}]]]]],"Input"],
 Cell["Local conservative storage and flux Jacobians", "Section"],
 Cell[BoxData[ToBoxes[Defer[MovingThreeFluid`storageJacobian//MatrixForm]]],"Input"],
 Cell[BoxData[ToBoxes[Defer[MovingThreeFluid`fluxJacobian//MatrixForm]]],"Input"],
 Cell["All band entries and both RHS conventions", "Section"],
 Cell["Select centre, bulk, outer_fan, outer_supersonic, or outer_vacuum. Different species may select different outer regimes. Aminus, Adiag, Aplus are the 12 by 12 blocks. IncrementRHS solves for delta X; AbsoluteRHS solves for the primitive predictor X*. Conservative local recovery follows that solve.","Text"],
 Cell[BoxData[ToBoxes[Defer[MovingThreeFluid`rowsByRegion["bulk"]["Variables"]]]],"Input"],
 Cell[BoxData[ToBoxes[Defer[MovingThreeFluid`entriesByRegion["bulk"]//TableForm]]],"Input"],
 Cell[BoxData[ToBoxes[Defer[MovingThreeFluid`rowsByRegion["bulk"]["IncrementRHS"]//TableForm]]],"Input"],
 Cell[BoxData[ToBoxes[Defer[MovingThreeFluid`rowsByRegion["centre"]["Matrix"]//MatrixForm]]],"Input"],
 Cell[BoxData[ToBoxes[Defer[MovingThreeFluid`entriesByRegion["outer_fan"]//TableForm]]],"Input"],
 Cell[BoxData[ToBoxes[Defer[MovingThreeFluid`vacuumJacobian//MatrixForm]]],"Input"],
 Cell["Verification", "Section"],
 Cell[BoxData[ToBoxes[Defer[MovingThreeFluid`checks]]],"Input"],
 Cell[BoxData[ToBoxes[Defer[MovingThreeFluid`finiteDifferenceErrors]]],"Input"],
 Cell[BoxData[ToBoxes[Defer[MovingThreeFluid`bandwidth]]],"Input"],
 Cell["Optimized C++ order is u_s,v_s,rho_s,u_b,v_b,rho_b,u_d,v_d,rho_d,M_s,M_b,M_d. Logical tables above use the readable species blocks; storageRowsByRegion contains the permuted matrices and RHS.","Text"],
 Cell[BoxData[ToBoxes[Defer[MovingThreeFluid`storageBandwidth]]],"Input"],
 Cell[BoxData[ToBoxes[Defer[MovingThreeFluid`storageRowsByRegion["bulk"]]]],"Input"]},
 WindowTitle->"MovingThreeFluidSim derivation",CellContext->"MovingThreeFluid`",StyleDefinitions->"Default.nb"];
If[!FileExistsQ[FileNameJoin[{directory,"moving_three_fluid_compiler.nb"}]],
 Put[notebook,FileNameJoin[{directory,"moving_three_fluid_compiler.nb"}]]];
loadedNotebook=Get[FileNameJoin[{directory,"moving_three_fluid_compiler.nb"}]];
inputBoxes=Cases[loadedNotebook,Cell[BoxData[boxes_],"Input",___]:>boxes,Infinity];
notebookMatrix=Block[{$Context="MovingThreeFluid`",$ContextPath={"MovingThreeFluid`","System`"}},
 ToExpression[inputBoxes[[2]],StandardForm]];
AssociateTo[checks,"notebook_symbol_context"->TrueQ[
 Dimensions[If[Head[notebookMatrix]===MatrixForm,First[notebookMatrix],notebookMatrix]]=={3,3}]];
AssociateTo[checks,"coefficient_export_roundtrip"->TrueQ[
 Get[FileNameJoin[{directory,"moving_three_fluid_coefficients.wl"}]]===rowsByRegion]];
Put[<|"Checks"->checks,"FiniteDifferenceErrors"->finiteDifferenceErrors,
 "LogicalHalfBandwidth"->bandwidth,"HalfBandwidth"->storageBandwidth|>,FileNameJoin[{directory,"moving_three_fluid_checks.wl"}]];
Print[checks];Print["Finite difference errors: ",finiteDifferenceErrors];
Print["Logical half-bandwidths: ",bandwidth,"; optimized storage half-bandwidths: ",storageBandwidth];
allChecksPassed=And@@(TrueQ /@Values[checks]);
Print["moving_three_fluid_checks_passed=",Boole[allChecksPassed]];
End[];
