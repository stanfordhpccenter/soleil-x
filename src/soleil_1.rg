struct particles_columns {
  cell : int3d;
  position : double[3];
  particle_velocity : double[3];
  density : double;
  particle_temperature : double;
  diameter : double;
  position_ghost : double[3];
  velocity_ghost : double[3];
  velocity_t_ghost : double[3];
  deltaVelocityOverRelaxationTime : double[3];
  deltaTemperatureTerm : double;
  position_old : double[3];
  velocity_old : double[3];
  temperature_old : double;
  position_new : double[3];
  velocity_new : double[3];
  temperature_new : double;
  position_t : double[3];
  velocity_t : double[3];
  temperature_t : double;
  __valid : bool;
}
struct Fluid_columns {
  rho : double;
  pressure : double;
  velocity : double[3];
  centerCoordinates : double[3];
  velocityGradientX : double[3];
  velocityGradientY : double[3];
  velocityGradientZ : double[3];
  temperature : double;
  rhoEnthalpy : double;
  kineticEnergy : double;
  sgsEnergy : double;
  sgsEddyViscosity : double;
  sgsEddyKappa : double;
  convectiveSpectralRadius : double;
  viscousSpectralRadius : double;
  heatConductionSpectralRadius : double;
  rhoVelocity : double[3];
  rhoEnergy : double;
  rhoBoundary : double;
  rhoVelocityBoundary : double[3];
  rhoEnergyBoundary : double;
  velocityBoundary : double[3];
  pressureBoundary : double;
  temperatureBoundary : double;
  velocityGradientXBoundary : double[3];
  velocityGradientYBoundary : double[3];
  velocityGradientZBoundary : double[3];
  rho_old : double;
  rhoVelocity_old : double[3];
  rhoEnergy_old : double;
  rho_new : double;
  rhoVelocity_new : double[3];
  rhoEnergy_new : double;
  rho_t : double;
  rhoVelocity_t : double[3];
  rhoEnergy_t : double;
  rhoFluxX : double;
  rhoVelocityFluxX : double[3];
  rhoEnergyFluxX : double;
  rhoFluxY : double;
  rhoVelocityFluxY : double[3];
  rhoEnergyFluxY : double;
  rhoFluxZ : double;
  rhoVelocityFluxZ : double[3];
  rhoEnergyFluxZ : double;
  PD : double;
  dissipation : double;
  dissipationFlux : double;
}
task InitParticlesUniform($particles#2 : region#1(ispace#1(int1d), particles_columns), $cells#4 : region#2(ispace#2(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($particles#2), writes($particles#2), reads($cells#4.velocity), reads($cells#4.centerCoordinates)
do
  var $pBase#6 : int32 = 0
  for $p#7 : int1d(particles_columns, $particles#2) in $particles#2 do
    $pBase#6 = int32($p#7)
    break
  end
  var $lo#8 : int3d = $cells#4.bounds.lo
  $lo#8.x = max($lo#8.x, 0)
  $lo#8.y = max($lo#8.y, 1)
  $lo#8.z = max($lo#8.z, 1)
  var $hi#9 : int3d = $cells#4.bounds.hi
  $hi#9.x = min($hi#9.x, ((378+0)-1))
  $hi#9.y = min($hi#9.y, ((831+1)-1))
  $hi#9.z = min($hi#9.z, ((732+1)-1))
  var $xSize#10 : int64 = (($hi#9.x-$lo#8.x)+1)
  var $ySize#11 : int64 = (($hi#9.y-$lo#8.y)+1)
  for $p#12 : int1d(particles_columns, $particles#2) in $particles#2 do
    if ((int32($p#12)-$pBase#6)<5) then
      (@$p#12).__valid = true
      var $relIdx#13 : int32 = (int32($p#12)-$pBase#6)
      var $c#14 : int3d = int3d({($lo#8.x+($relIdx#13%$xSize#10)), ($lo#8.y+(($relIdx#13/$xSize#10)%$ySize#11)), ($lo#8.z+(($relIdx#13/$xSize#10)/$ySize#11))})
      (@$p#12).cell = $c#14
      (@$p#12).position = $cells#4[(@$p#12).cell].centerCoordinates
      (@$p#12).particle_velocity = $cells#4[(@$p#12).cell].velocity
      (@$p#12).density = 0.4769553
      (@$p#12).particle_temperature = 0.7490512
      (@$p#12).diameter = 0.5250806
    else
    end
  end
end
task AddRadiation($particles#191 : region#92(ispace#92(int1d), particles_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($particles#191.density), reads($particles#191.diameter), reads($particles#191.temperature_t), writes($particles#191.temperature_t)
do
  for $p#193 : int1d(particles_columns, $particles#191) in $particles#191 do
    var $crossSectionArea#194 : double = (((2*acos(0))*pow((@$p#193).diameter, 2))/4)
    var $volume#195 : double = (((2*acos(0))*pow((@$p#193).diameter, 3))/6)
    var $mass#196 : double = ($volume#195*(@$p#193).density)
    var $absorbedRadiationIntensity#197 : double = ((0.6796266*0.258018)*$crossSectionArea#194)
    (@$p#193).temperature_t += ($absorbedRadiationIntensity#197/($mass#196*0.3767944))
  end
end
task particles_initValidField($r#1127 : region#139(ispace#139(int1d), particles_columns))
-- leaf (false), inner (false), idempotent (false)
where
  writes($r#1127.__valid)
do
  for $e#1129 : int1d(particles_columns, $r#1127) in $r#1127 do
    (@$e#1129).__valid = false
  end
end
task Flow_InitializeCell($dom#1156 : region#152(ispace#152(int3d), Fluid_columns), $Fluid#1158 : region#153(ispace#153(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#1158.PD), writes($Fluid#1158.PD), reads($Fluid#1158.centerCoordinates), writes($Fluid#1158.centerCoordinates), reads($Fluid#1158.convectiveSpectralRadius), writes($Fluid#1158.convectiveSpectralRadius), reads($Fluid#1158.dissipation), writes($Fluid#1158.dissipation), reads($Fluid#1158.dissipationFlux), writes($Fluid#1158.dissipationFlux), reads($Fluid#1158.heatConductionSpectralRadius), writes($Fluid#1158.heatConductionSpectralRadius), reads($Fluid#1158.kineticEnergy), writes($Fluid#1158.kineticEnergy), reads($Fluid#1158.pressure), writes($Fluid#1158.pressure), reads($Fluid#1158.pressureBoundary), writes($Fluid#1158.pressureBoundary), reads($Fluid#1158.rho), writes($Fluid#1158.rho), reads($Fluid#1158.rhoBoundary), writes($Fluid#1158.rhoBoundary), reads($Fluid#1158.rhoEnergy), writes($Fluid#1158.rhoEnergy), reads($Fluid#1158.rhoEnergyBoundary), writes($Fluid#1158.rhoEnergyBoundary), reads($Fluid#1158.rhoEnergyFluxX), writes($Fluid#1158.rhoEnergyFluxX), reads($Fluid#1158.rhoEnergyFluxY), writes($Fluid#1158.rhoEnergyFluxY), reads($Fluid#1158.rhoEnergyFluxZ), writes($Fluid#1158.rhoEnergyFluxZ), reads($Fluid#1158.rhoEnergy_new), writes($Fluid#1158.rhoEnergy_new), reads($Fluid#1158.rhoEnergy_old), writes($Fluid#1158.rhoEnergy_old), reads($Fluid#1158.rhoEnergy_t), writes($Fluid#1158.rhoEnergy_t), reads($Fluid#1158.rhoEnthalpy), writes($Fluid#1158.rhoEnthalpy), reads($Fluid#1158.rhoFluxX), writes($Fluid#1158.rhoFluxX), reads($Fluid#1158.rhoFluxY), writes($Fluid#1158.rhoFluxY), reads($Fluid#1158.rhoFluxZ), writes($Fluid#1158.rhoFluxZ), reads($Fluid#1158.rhoVelocity), writes($Fluid#1158.rhoVelocity), reads($Fluid#1158.rhoVelocityBoundary), writes($Fluid#1158.rhoVelocityBoundary), reads($Fluid#1158.rhoVelocityFluxX), writes($Fluid#1158.rhoVelocityFluxX), reads($Fluid#1158.rhoVelocityFluxY), writes($Fluid#1158.rhoVelocityFluxY), reads($Fluid#1158.rhoVelocityFluxZ), writes($Fluid#1158.rhoVelocityFluxZ), reads($Fluid#1158.rhoVelocity_new), writes($Fluid#1158.rhoVelocity_new), reads($Fluid#1158.rhoVelocity_old), writes($Fluid#1158.rhoVelocity_old), reads($Fluid#1158.rhoVelocity_t), writes($Fluid#1158.rhoVelocity_t), reads($Fluid#1158.rho_new), writes($Fluid#1158.rho_new), reads($Fluid#1158.rho_old), writes($Fluid#1158.rho_old), reads($Fluid#1158.rho_t), writes($Fluid#1158.rho_t), reads($Fluid#1158.sgsEddyKappa), writes($Fluid#1158.sgsEddyKappa), reads($Fluid#1158.sgsEddyViscosity), writes($Fluid#1158.sgsEddyViscosity), reads($Fluid#1158.sgsEnergy), writes($Fluid#1158.sgsEnergy), reads($Fluid#1158.temperature), writes($Fluid#1158.temperature), reads($Fluid#1158.temperatureBoundary), writes($Fluid#1158.temperatureBoundary), reads($Fluid#1158.velocity), writes($Fluid#1158.velocity), reads($Fluid#1158.velocityBoundary), writes($Fluid#1158.velocityBoundary), reads($Fluid#1158.velocityGradientX), writes($Fluid#1158.velocityGradientX), reads($Fluid#1158.velocityGradientXBoundary), writes($Fluid#1158.velocityGradientXBoundary), reads($Fluid#1158.velocityGradientY), writes($Fluid#1158.velocityGradientY), reads($Fluid#1158.velocityGradientYBoundary), writes($Fluid#1158.velocityGradientYBoundary), reads($Fluid#1158.velocityGradientZ), writes($Fluid#1158.velocityGradientZ), reads($Fluid#1158.velocityGradientZBoundary), writes($Fluid#1158.velocityGradientZBoundary), reads($Fluid#1158.viscousSpectralRadius), writes($Fluid#1158.viscousSpectralRadius), $dom#1156 <= $Fluid#1158
do
  for $c#1161 : int3d(Fluid_columns, $dom#1156) in $dom#1156 do
    $Fluid#1158[$c#1161].rho = double(0)
    $Fluid#1158[$c#1161].pressure = double(0)
    $Fluid#1158[$c#1161].velocity = double[3](array(double(0), double(0), double(0)))
    $Fluid#1158[$c#1161].centerCoordinates = double[3](array(double(0), double(0), double(0)))
    $Fluid#1158[$c#1161].velocityGradientX = double[3](array(double(0), double(0), double(0)))
    $Fluid#1158[$c#1161].velocityGradientY = double[3](array(double(0), double(0), double(0)))
    $Fluid#1158[$c#1161].velocityGradientZ = double[3](array(double(0), double(0), double(0)))
    $Fluid#1158[$c#1161].temperature = double(0)
    $Fluid#1158[$c#1161].rhoEnthalpy = double(0)
    $Fluid#1158[$c#1161].kineticEnergy = double(0)
    $Fluid#1158[$c#1161].sgsEnergy = double(0)
    $Fluid#1158[$c#1161].sgsEddyViscosity = double(0)
    $Fluid#1158[$c#1161].sgsEddyKappa = double(0)
    $Fluid#1158[$c#1161].convectiveSpectralRadius = double(0)
    $Fluid#1158[$c#1161].viscousSpectralRadius = double(0)
    $Fluid#1158[$c#1161].heatConductionSpectralRadius = double(0)
    $Fluid#1158[$c#1161].rhoVelocity = double[3](array(double(0), double(0), double(0)))
    $Fluid#1158[$c#1161].rhoEnergy = double(0)
    $Fluid#1158[$c#1161].rhoBoundary = double(0)
    $Fluid#1158[$c#1161].rhoVelocityBoundary = double[3](array(double(0), double(0), double(0)))
    $Fluid#1158[$c#1161].rhoEnergyBoundary = double(0)
    $Fluid#1158[$c#1161].velocityBoundary = double[3](array(double(0), double(0), double(0)))
    $Fluid#1158[$c#1161].pressureBoundary = double(0)
    $Fluid#1158[$c#1161].temperatureBoundary = double(0)
    $Fluid#1158[$c#1161].velocityGradientXBoundary = double[3](array(double(0), double(0), double(0)))
    $Fluid#1158[$c#1161].velocityGradientYBoundary = double[3](array(double(0), double(0), double(0)))
    $Fluid#1158[$c#1161].velocityGradientZBoundary = double[3](array(double(0), double(0), double(0)))
    $Fluid#1158[$c#1161].rho_old = double(0)
    $Fluid#1158[$c#1161].rhoVelocity_old = double[3](array(double(0), double(0), double(0)))
    $Fluid#1158[$c#1161].rhoEnergy_old = double(0)
    $Fluid#1158[$c#1161].rho_new = double(0)
    $Fluid#1158[$c#1161].rhoVelocity_new = double[3](array(double(0), double(0), double(0)))
    $Fluid#1158[$c#1161].rhoEnergy_new = double(0)
    $Fluid#1158[$c#1161].rho_t = double(0)
    $Fluid#1158[$c#1161].rhoVelocity_t = double[3](array(double(0), double(0), double(0)))
    $Fluid#1158[$c#1161].rhoEnergy_t = double(0)
    $Fluid#1158[$c#1161].rhoFluxX = double(0)
    $Fluid#1158[$c#1161].rhoVelocityFluxX = double[3](array(double(0), double(0), double(0)))
    $Fluid#1158[$c#1161].rhoEnergyFluxX = double(0)
    $Fluid#1158[$c#1161].rhoFluxY = double(0)
    $Fluid#1158[$c#1161].rhoVelocityFluxY = double[3](array(double(0), double(0), double(0)))
    $Fluid#1158[$c#1161].rhoEnergyFluxY = double(0)
    $Fluid#1158[$c#1161].rhoFluxZ = double(0)
    $Fluid#1158[$c#1161].rhoVelocityFluxZ = double[3](array(double(0), double(0), double(0)))
    $Fluid#1158[$c#1161].rhoEnergyFluxZ = double(0)
    $Fluid#1158[$c#1161].PD = double(0)
    $Fluid#1158[$c#1161].dissipation = double(0)
    $Fluid#1158[$c#1161].dissipationFlux = double(0)
  end
end
task Flow_InitializeCenterCoordinates($dom#1287 : region#163(ispace#163(int3d), Fluid_columns), $Fluid#1289 : region#164(ispace#164(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#1289.centerCoordinates), writes($Fluid#1289.centerCoordinates), $dom#1287 <= $Fluid#1289
do
  for $c#1296 : int3d(Fluid_columns, $dom#1287) in $dom#1287 do
    var $xy#1297 : double[3] = double[3](array((double(0.6987052)+((double(0.5988125)/double(int32(378)))*(double(int3d($c#1296).x)+double(0.5)))), (double(0.99700342743682)+((double(0.84872084512635)/double(int32(833)))*(double(int3d($c#1296).y)+double(0.5)))), (double(0.70456500505464)+((double(0.83330648989071)/double(int32(734)))*(double(int3d($c#1296).z)+double(0.5))))))
    $Fluid#1289[$c#1296].centerCoordinates = double[3](array(double($xy#1297[int32(0)]), double($xy#1297[int32(1)]), double($xy#1297[int32(2)])))
  end
end
task Flow_InitializeUniform($dom#1336 : region#176(ispace#176(int3d), Fluid_columns), $Fluid#1338 : region#177(ispace#177(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#1338.pressure), writes($Fluid#1338.pressure), reads($Fluid#1338.rho), writes($Fluid#1338.rho), reads($Fluid#1338.velocity), writes($Fluid#1338.velocity), $dom#1336 <= $Fluid#1338
do
  for $c#1341 : int3d(Fluid_columns, $dom#1336) in $dom#1336 do
    $Fluid#1338[$c#1341].rho = array(double(0.4472728), double(0.4001904), double(0.3713178), double(0.6807681), double(0.4900506))[int32(0)]
    $Fluid#1338[$c#1341].pressure = array(double(0.4472728), double(0.4001904), double(0.3713178), double(0.6807681), double(0.4900506))[int32(1)]
    $Fluid#1338[$c#1341].velocity[int32(0)] = array(double(0.4472728), double(0.4001904), double(0.3713178), double(0.6807681), double(0.4900506))[int32(2)]
    $Fluid#1338[$c#1341].velocity[int32(1)] = array(double(0.4472728), double(0.4001904), double(0.3713178), double(0.6807681), double(0.4900506))[int32(3)]
    $Fluid#1338[$c#1341].velocity[int32(2)] = array(double(0.4472728), double(0.4001904), double(0.3713178), double(0.6807681), double(0.4900506))[int32(4)]
  end
end
terra vs_mul_double_3($a : double[3],$b : double) : double[3]
    return array([&double]($a)[0] * $b, [&double]($a)[1] * $b, [&double]($a)[2] * $b)
end
terra dot_double_3($a : double[3],$b : double[3]) : double
    return [&double]($a)[0] * [&double]($b)[0] + [&double]($a)[1] * [&double]($b)[1] + [&double]($a)[2] * [&double]($b)[2]
end
task Flow_UpdateConservedFromPrimitive($dom#1415 : region#187(ispace#187(int3d), Fluid_columns), $Fluid#1417 : region#188(ispace#188(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#1417.pressure), reads($Fluid#1417.rho), reads($Fluid#1417.rhoEnergy), writes($Fluid#1417.rhoEnergy), reads($Fluid#1417.rhoVelocity), writes($Fluid#1417.rhoVelocity), reads($Fluid#1417.sgsEnergy), reads($Fluid#1417.velocity), $dom#1415 <= $Fluid#1417
do
  for $c#1438 : int3d(Fluid_columns, $dom#1415) in $dom#1415 do
    if (not ((((((max(int32((uint64(int32(0))-int3d($c#1438).x)), int32(0))>int32(0)) or (max(int32((int3d($c#1438).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#1438).y)), int32(0))>int32(0))) or (max(int32((int3d($c#1438).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#1438).z)), int32(0))>int32(0))) or (max(int32((int3d($c#1438).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0)))) then
      var $tmpTemperature#1439 : double = ($Fluid#1417[$c#1438].pressure/(double(0.7272462)*$Fluid#1417[$c#1438].rho))
      var $velocity#1440 : double[3] = $Fluid#1417[$c#1438].velocity
      $Fluid#1417[$c#1438].rhoVelocity = vs_mul_double_3($Fluid#1417[$c#1438].velocity, $Fluid#1417[$c#1438].rho)
      var $cv#1441 : double = (double(0.7272462)/(double(0.8198141)-double(1)))
      $Fluid#1417[$c#1438].rhoEnergy = (($Fluid#1417[$c#1438].rho*(($cv#1441*$tmpTemperature#1439)+(double(0.5)*dot_double_3($velocity#1440, $velocity#1440))))+$Fluid#1417[$c#1438].sgsEnergy)
    else
    end
  end
end
terra vs_div_double_3($a : double[3],$b : double) : double[3]
    return array([&double]($a)[0] / $b, [&double]($a)[1] / $b, [&double]($a)[2] / $b)
end
task Flow_UpdateAuxiliaryVelocity($dom#1480 : region#204(ispace#204(int3d), Fluid_columns), $Fluid#1482 : region#205(ispace#205(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#1482.kineticEnergy), writes($Fluid#1482.kineticEnergy), reads($Fluid#1482.rho), reads($Fluid#1482.rhoVelocity), reads($Fluid#1482.velocity), writes($Fluid#1482.velocity), $dom#1480 <= $Fluid#1482
do
  for $c#1491 : int3d(Fluid_columns, $dom#1480) in $dom#1480 do
    if (not ((((((max(int32((uint64(int32(0))-int3d($c#1491).x)), int32(0))>int32(0)) or (max(int32((int3d($c#1491).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#1491).y)), int32(0))>int32(0))) or (max(int32((int3d($c#1491).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#1491).z)), int32(0))>int32(0))) or (max(int32((int3d($c#1491).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0)))) then
      var $velocity#1492 : double[3] = vs_div_double_3($Fluid#1482[$c#1491].rhoVelocity, $Fluid#1482[$c#1491].rho)
      $Fluid#1482[$c#1491].velocity = $velocity#1492
      $Fluid#1482[$c#1491].kineticEnergy = ((double(0.5)*$Fluid#1482[$c#1491].rho)*dot_double_3($velocity#1492, $velocity#1492))
    else
    end
  end
end
terra vv_mul_double_3($a : double[3],$b : double[3]) : double[3]
    return array([&double]($a)[0] * [&double]($b)[0], [&double]($a)[1] * [&double]($b)[1], [&double]($a)[2] * [&double]($b)[2])
end
terra vv_add_double_3($a : double[3],$b : double[3]) : double[3]
    return array([&double]($a)[0] + [&double]($b)[0], [&double]($a)[1] + [&double]($b)[1], [&double]($a)[2] + [&double]($b)[2])
end
task Flow_UpdateGhostConservedStep1($dom#1525 : region#217(ispace#217(int3d), Fluid_columns), $Fluid#1527 : region#218(ispace#218(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#1527.pressure), reads($Fluid#1527.rho), reads($Fluid#1527.rhoBoundary), writes($Fluid#1527.rhoBoundary), reads($Fluid#1527.rhoEnergyBoundary), writes($Fluid#1527.rhoEnergyBoundary), reads($Fluid#1527.rhoVelocity), reads($Fluid#1527.rhoVelocityBoundary), writes($Fluid#1527.rhoVelocityBoundary), reads($Fluid#1527.temperature), $dom#1525 <= $Fluid#1527
do
  for $c#1992 : int3d(Fluid_columns, $dom#1525) in $dom#1525 do
    if (max(int32((uint64(int32(0))-int3d($c#1992).x)), int32(0))>int32(0)) then
      var $c_bnd#1993 : int3d(Fluid_columns, $dom#1525) = $c#1992
      var $c_int#1994 : int3d = (($c#1992+{1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#1995 : double[3] = array(double(1), double(1), double(1))
      var $bnd_velocity#1996 : double[3] = array(double(0), double(0), double(0))
      var $bnd_temperature#1997 : double = double(-1)
      var $rho#1998 : double = double(double(0))
      var $temp_wall#1999 : double = double(double(0))
      var $temperature#2000 : double = double(double(0))
      var $velocity#2001 : double[3] = double[3](array(double(0), double(0), double(0)))
      var $cv#2002 : double = (double(0.7272462)/(double(0.8198141)-double(1)))
      var $velocity#2003 : double[3] = double[3](array(double(0), double(0), double(0)))
      $velocity#2003 = vv_add_double_3(vv_mul_double_3(vs_div_double_3($Fluid#1527[$c_int#1994].rhoVelocity, $Fluid#1527[$c_int#1994].rho), $sign#1995), $bnd_velocity#1996)
      $temp_wall#1999 = $Fluid#1527[$c_int#1994].temperature
      if ($bnd_temperature#1997>double(0)) then
        $temp_wall#1999 = $bnd_temperature#1997
      else
      end
      $temperature#2000 = ((double(2)*$temp_wall#1999)-$Fluid#1527[$c_int#1994].temperature)
      $rho#1998 = ($Fluid#1527[$c_int#1994].pressure/(double(0.7272462)*$temperature#2000))
      $Fluid#1527[$c_bnd#1993].rhoBoundary = $rho#1998
      $Fluid#1527[$c_bnd#1993].rhoVelocityBoundary = vs_mul_double_3($velocity#2003, $rho#1998)
      $Fluid#1527[$c_bnd#1993].rhoEnergyBoundary = ($rho#1998*(($cv#2002*$temperature#2000)+(double(0.5)*dot_double_3($velocity#2003, $velocity#2003))))
    else
    end
    if (max(int32((int3d($c#1992).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0)) then
      var $c_bnd#2004 : int3d(Fluid_columns, $dom#1525) = $c#1992
      var $c_int#2005 : int3d = (($c#1992+{-1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#2006 : double[3] = array(double(1), double(1), double(1))
      var $bnd_velocity#2007 : double[3] = array(double(0), double(0), double(0))
      var $bnd_temperature#2008 : double = double(-1)
      var $rho#2009 : double = double(double(0))
      var $temp_wall#2010 : double = double(double(0))
      var $temperature#2011 : double = double(double(0))
      var $velocity#2012 : double[3] = double[3](array(double(0), double(0), double(0)))
      var $cv#2013 : double = (double(0.7272462)/(double(0.8198141)-double(1)))
      var $velocity#2014 : double[3] = double[3](array(double(0), double(0), double(0)))
      $velocity#2014 = vv_add_double_3(vv_mul_double_3(vs_div_double_3($Fluid#1527[$c_int#2005].rhoVelocity, $Fluid#1527[$c_int#2005].rho), $sign#2006), $bnd_velocity#2007)
      $temp_wall#2010 = $Fluid#1527[$c_int#2005].temperature
      if ($bnd_temperature#2008>double(0)) then
        $temp_wall#2010 = $bnd_temperature#2008
      else
      end
      $temperature#2011 = ((double(2)*$temp_wall#2010)-$Fluid#1527[$c_int#2005].temperature)
      $rho#2009 = ($Fluid#1527[$c_int#2005].pressure/(double(0.7272462)*$temperature#2011))
      $Fluid#1527[$c_bnd#2004].rhoBoundary = $rho#2009
      $Fluid#1527[$c_bnd#2004].rhoVelocityBoundary = vs_mul_double_3($velocity#2014, $rho#2009)
      $Fluid#1527[$c_bnd#2004].rhoEnergyBoundary = ($rho#2009*(($cv#2013*$temperature#2011)+(double(0.5)*dot_double_3($velocity#2014, $velocity#2014))))
    else
    end
    if (max(int32((uint64(int32(1))-int3d($c#1992).y)), int32(0))>int32(0)) then
      var $c_bnd#2015 : int3d(Fluid_columns, $dom#1525) = $c#1992
      var $c_int#2016 : int3d = (($c#1992+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#2017 : double[3] = array(double(-1), double(-1), double(-1))
      var $bnd_velocity#2018 : double[3] = array(double(1.9836594), double(1.7248182), double(0.6582188))
      var $bnd_temperature#2019 : double = double(-1)
      var $rho#2020 : double = double(double(0))
      var $temp_wall#2021 : double = double(double(0))
      var $temperature#2022 : double = double(double(0))
      var $velocity#2023 : double[3] = double[3](array(double(0), double(0), double(0)))
      var $cv#2024 : double = (double(0.7272462)/(double(0.8198141)-double(1)))
      var $velocity#2025 : double[3] = double[3](array(double(0), double(0), double(0)))
      $velocity#2025 = vv_add_double_3(vv_mul_double_3(vs_div_double_3($Fluid#1527[$c_int#2016].rhoVelocity, $Fluid#1527[$c_int#2016].rho), $sign#2017), $bnd_velocity#2018)
      $temp_wall#2021 = $Fluid#1527[$c_int#2016].temperature
      if ($bnd_temperature#2019>double(0)) then
        $temp_wall#2021 = $bnd_temperature#2019
      else
      end
      $temperature#2022 = ((double(2)*$temp_wall#2021)-$Fluid#1527[$c_int#2016].temperature)
      $rho#2020 = ($Fluid#1527[$c_int#2016].pressure/(double(0.7272462)*$temperature#2022))
      $Fluid#1527[$c_bnd#2015].rhoBoundary = $rho#2020
      $Fluid#1527[$c_bnd#2015].rhoVelocityBoundary = vs_mul_double_3($velocity#2025, $rho#2020)
      $Fluid#1527[$c_bnd#2015].rhoEnergyBoundary = ($rho#2020*(($cv#2024*$temperature#2022)+(double(0.5)*dot_double_3($velocity#2025, $velocity#2025))))
    else
    end
    if (max(int32((int3d($c#1992).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0)) then
      var $c_bnd#2026 : int3d(Fluid_columns, $dom#1525) = $c#1992
      var $c_int#2027 : int3d = (($c#1992+{0, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#2028 : double[3] = array(double(-1), double(-1), double(-1))
      var $bnd_velocity#2029 : double[3] = array(double(0.8413202), double(1.3283396), double(1.134841))
      var $bnd_temperature#2030 : double = double(-1)
      var $rho#2031 : double = double(double(0))
      var $temp_wall#2032 : double = double(double(0))
      var $temperature#2033 : double = double(double(0))
      var $velocity#2034 : double[3] = double[3](array(double(0), double(0), double(0)))
      var $cv#2035 : double = (double(0.7272462)/(double(0.8198141)-double(1)))
      var $velocity#2036 : double[3] = double[3](array(double(0), double(0), double(0)))
      $velocity#2036 = vv_add_double_3(vv_mul_double_3(vs_div_double_3($Fluid#1527[$c_int#2027].rhoVelocity, $Fluid#1527[$c_int#2027].rho), $sign#2028), $bnd_velocity#2029)
      $temp_wall#2032 = $Fluid#1527[$c_int#2027].temperature
      if ($bnd_temperature#2030>double(0)) then
        $temp_wall#2032 = $bnd_temperature#2030
      else
      end
      $temperature#2033 = ((double(2)*$temp_wall#2032)-$Fluid#1527[$c_int#2027].temperature)
      $rho#2031 = ($Fluid#1527[$c_int#2027].pressure/(double(0.7272462)*$temperature#2033))
      $Fluid#1527[$c_bnd#2026].rhoBoundary = $rho#2031
      $Fluid#1527[$c_bnd#2026].rhoVelocityBoundary = vs_mul_double_3($velocity#2036, $rho#2031)
      $Fluid#1527[$c_bnd#2026].rhoEnergyBoundary = ($rho#2031*(($cv#2035*$temperature#2033)+(double(0.5)*dot_double_3($velocity#2036, $velocity#2036))))
    else
    end
    if (max(int32((uint64(int32(1))-int3d($c#1992).z)), int32(0))>int32(0)) then
      var $c_bnd#2037 : int3d(Fluid_columns, $dom#1525) = $c#1992
      var $c_int#2038 : int3d = (($c#1992+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#2039 : double[3] = array(double(1), double(1), double(-1))
      var $bnd_velocity#2040 : double[3] = array(double(0), double(0), double(0))
      var $bnd_temperature#2041 : double = double(-1)
      var $rho#2042 : double = double(double(0))
      var $temp_wall#2043 : double = double(double(0))
      var $temperature#2044 : double = double(double(0))
      var $velocity#2045 : double[3] = double[3](array(double(0), double(0), double(0)))
      var $cv#2046 : double = (double(0.7272462)/(double(0.8198141)-double(1)))
      var $velocity#2047 : double[3] = double[3](array(double(0), double(0), double(0)))
      $velocity#2047 = vv_add_double_3(vv_mul_double_3(vs_div_double_3($Fluid#1527[$c_int#2038].rhoVelocity, $Fluid#1527[$c_int#2038].rho), $sign#2039), $bnd_velocity#2040)
      $temp_wall#2043 = $Fluid#1527[$c_int#2038].temperature
      if ($bnd_temperature#2041>double(0)) then
        $temp_wall#2043 = $bnd_temperature#2041
      else
      end
      $temperature#2044 = ((double(2)*$temp_wall#2043)-$Fluid#1527[$c_int#2038].temperature)
      $rho#2042 = ($Fluid#1527[$c_int#2038].pressure/(double(0.7272462)*$temperature#2044))
      $Fluid#1527[$c_bnd#2037].rhoBoundary = $rho#2042
      $Fluid#1527[$c_bnd#2037].rhoVelocityBoundary = vs_mul_double_3($velocity#2047, $rho#2042)
      $Fluid#1527[$c_bnd#2037].rhoEnergyBoundary = ($rho#2042*(($cv#2046*$temperature#2044)+(double(0.5)*dot_double_3($velocity#2047, $velocity#2047))))
    else
    end
    if (max(int32((int3d($c#1992).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0)) then
      var $c_bnd#2048 : int3d(Fluid_columns, $dom#1525) = $c#1992
      var $c_int#2049 : int3d = (($c#1992+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#2050 : double[3] = array(double(1), double(1), double(-1))
      var $bnd_velocity#2051 : double[3] = array(double(0), double(0), double(0))
      var $bnd_temperature#2052 : double = double(-1)
      var $rho#2053 : double = double(double(0))
      var $temp_wall#2054 : double = double(double(0))
      var $temperature#2055 : double = double(double(0))
      var $velocity#2056 : double[3] = double[3](array(double(0), double(0), double(0)))
      var $cv#2057 : double = (double(0.7272462)/(double(0.8198141)-double(1)))
      var $velocity#2058 : double[3] = double[3](array(double(0), double(0), double(0)))
      $velocity#2058 = vv_add_double_3(vv_mul_double_3(vs_div_double_3($Fluid#1527[$c_int#2049].rhoVelocity, $Fluid#1527[$c_int#2049].rho), $sign#2050), $bnd_velocity#2051)
      $temp_wall#2054 = $Fluid#1527[$c_int#2049].temperature
      if ($bnd_temperature#2052>double(0)) then
        $temp_wall#2054 = $bnd_temperature#2052
      else
      end
      $temperature#2055 = ((double(2)*$temp_wall#2054)-$Fluid#1527[$c_int#2049].temperature)
      $rho#2053 = ($Fluid#1527[$c_int#2049].pressure/(double(0.7272462)*$temperature#2055))
      $Fluid#1527[$c_bnd#2048].rhoBoundary = $rho#2053
      $Fluid#1527[$c_bnd#2048].rhoVelocityBoundary = vs_mul_double_3($velocity#2058, $rho#2053)
      $Fluid#1527[$c_bnd#2048].rhoEnergyBoundary = ($rho#2053*(($cv#2057*$temperature#2055)+(double(0.5)*dot_double_3($velocity#2058, $velocity#2058))))
    else
    end
  end
end
task Flow_UpdateGhostConservedStep2($dom#2658 : region#432(ispace#432(int3d), Fluid_columns), $Fluid#2660 : region#433(ispace#433(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#2660.rho), writes($Fluid#2660.rho), reads($Fluid#2660.rhoBoundary), reads($Fluid#2660.rhoEnergy), writes($Fluid#2660.rhoEnergy), reads($Fluid#2660.rhoEnergyBoundary), reads($Fluid#2660.rhoVelocity), writes($Fluid#2660.rhoVelocity), reads($Fluid#2660.rhoVelocityBoundary), $dom#2658 <= $Fluid#2660
do
  for $c#2663 : int3d(Fluid_columns, $dom#2658) in $dom#2658 do
    if ((((((max(int32((uint64(int32(0))-int3d($c#2663).x)), int32(0))>int32(0)) or (max(int32((int3d($c#2663).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#2663).y)), int32(0))>int32(0))) or (max(int32((int3d($c#2663).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#2663).z)), int32(0))>int32(0))) or (max(int32((int3d($c#2663).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0))) then
      $Fluid#2660[$c#2663].rho = $Fluid#2660[$c#2663].rhoBoundary
      $Fluid#2660[$c#2663].rhoVelocity = $Fluid#2660[$c#2663].rhoVelocityBoundary
      $Fluid#2660[$c#2663].rhoEnergy = $Fluid#2660[$c#2663].rhoEnergyBoundary
    else
    end
  end
end
task Flow_UpdateGhostVelocityStep1($dom#2687 : region#443(ispace#443(int3d), Fluid_columns), $Fluid#2689 : region#444(ispace#444(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#2689.velocity), reads($Fluid#2689.velocityBoundary), writes($Fluid#2689.velocityBoundary), $dom#2687 <= $Fluid#2689
do
  for $c#2860 : int3d(Fluid_columns, $dom#2687) in $dom#2687 do
    if (max(int32((uint64(int32(0))-int3d($c#2860).x)), int32(0))>int32(0)) then
      var $c_bnd#2861 : int3d(Fluid_columns, $dom#2687) = $c#2860
      var $c_int#2862 : int3d = (($c#2860+{1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#2863 : double[3] = array(double(1), double(1), double(1))
      var $bnd_velocity#2864 : double[3] = array(double(0), double(0), double(0))
      $Fluid#2689[$c_bnd#2861].velocityBoundary = vv_add_double_3(vv_mul_double_3($Fluid#2689[$c_int#2862].velocity, $sign#2863), $bnd_velocity#2864)
    else
    end
    if (max(int32((int3d($c#2860).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0)) then
      var $c_bnd#2865 : int3d(Fluid_columns, $dom#2687) = $c#2860
      var $c_int#2866 : int3d = (($c#2860+{-1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#2867 : double[3] = array(double(1), double(1), double(1))
      var $bnd_velocity#2868 : double[3] = array(double(0), double(0), double(0))
      $Fluid#2689[$c_bnd#2865].velocityBoundary = vv_add_double_3(vv_mul_double_3($Fluid#2689[$c_int#2866].velocity, $sign#2867), $bnd_velocity#2868)
    else
    end
    if (max(int32((uint64(int32(1))-int3d($c#2860).y)), int32(0))>int32(0)) then
      var $c_bnd#2869 : int3d(Fluid_columns, $dom#2687) = $c#2860
      var $c_int#2870 : int3d = (($c#2860+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#2871 : double[3] = array(double(-1), double(-1), double(-1))
      var $bnd_velocity#2872 : double[3] = array(double(1.9836594), double(1.7248182), double(0.6582188))
      $Fluid#2689[$c_bnd#2869].velocityBoundary = vv_add_double_3(vv_mul_double_3($Fluid#2689[$c_int#2870].velocity, $sign#2871), $bnd_velocity#2872)
    else
    end
    if (max(int32((int3d($c#2860).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0)) then
      var $c_bnd#2873 : int3d(Fluid_columns, $dom#2687) = $c#2860
      var $c_int#2874 : int3d = (($c#2860+{0, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#2875 : double[3] = array(double(-1), double(-1), double(-1))
      var $bnd_velocity#2876 : double[3] = array(double(0.8413202), double(1.3283396), double(1.134841))
      $Fluid#2689[$c_bnd#2873].velocityBoundary = vv_add_double_3(vv_mul_double_3($Fluid#2689[$c_int#2874].velocity, $sign#2875), $bnd_velocity#2876)
    else
    end
    if (max(int32((uint64(int32(1))-int3d($c#2860).z)), int32(0))>int32(0)) then
      var $c_bnd#2877 : int3d(Fluid_columns, $dom#2687) = $c#2860
      var $c_int#2878 : int3d = (($c#2860+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#2879 : double[3] = array(double(1), double(1), double(-1))
      var $bnd_velocity#2880 : double[3] = array(double(0), double(0), double(0))
      $Fluid#2689[$c_bnd#2877].velocityBoundary = vv_add_double_3(vv_mul_double_3($Fluid#2689[$c_int#2878].velocity, $sign#2879), $bnd_velocity#2880)
    else
    end
    if (max(int32((int3d($c#2860).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0)) then
      var $c_bnd#2881 : int3d(Fluid_columns, $dom#2687) = $c#2860
      var $c_int#2882 : int3d = (($c#2860+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#2883 : double[3] = array(double(1), double(1), double(-1))
      var $bnd_velocity#2884 : double[3] = array(double(0), double(0), double(0))
      $Fluid#2689[$c_bnd#2881].velocityBoundary = vv_add_double_3(vv_mul_double_3($Fluid#2689[$c_int#2882].velocity, $sign#2883), $bnd_velocity#2884)
    else
    end
  end
end
task Flow_UpdateGhostVelocityStep2($dom#3142 : region#526(ispace#526(int3d), Fluid_columns), $Fluid#3144 : region#527(ispace#527(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#3144.velocity), writes($Fluid#3144.velocity), reads($Fluid#3144.velocityBoundary), $dom#3142 <= $Fluid#3144
do
  for $c#3147 : int3d(Fluid_columns, $dom#3142) in $dom#3142 do
    if ((((((max(int32((uint64(int32(0))-int3d($c#3147).x)), int32(0))>int32(0)) or (max(int32((int3d($c#3147).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#3147).y)), int32(0))>int32(0))) or (max(int32((int3d($c#3147).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#3147).z)), int32(0))>int32(0))) or (max(int32((int3d($c#3147).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0))) then
      $Fluid#3144[$c#3147].velocity = $Fluid#3144[$c#3147].velocityBoundary
    else
    end
  end
end
terra vv_sub_double_3($a : double[3],$b : double[3]) : double[3]
    return array([&double]($a)[0] - [&double]($b)[0], [&double]($a)[1] - [&double]($b)[1], [&double]($a)[2] - [&double]($b)[2])
end
task Flow_ComputeVelocityGradientAll($dom#3171 : region#537(ispace#537(int3d), Fluid_columns), $Fluid#3173 : region#538(ispace#538(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#3173.velocity), reads($Fluid#3173.velocityGradientX), writes($Fluid#3173.velocityGradientX), reads($Fluid#3173.velocityGradientY), writes($Fluid#3173.velocityGradientY), reads($Fluid#3173.velocityGradientZ), writes($Fluid#3173.velocityGradientZ), $dom#3171 <= $Fluid#3173
do
  for $c#3176 : int3d(Fluid_columns, $dom#3171) in $dom#3171 do
    if (not ((((((max(int32((uint64(int32(0))-int3d($c#3176).x)), int32(0))>int32(0)) or (max(int32((int3d($c#3176).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#3176).y)), int32(0))>int32(0))) or (max(int32((int3d($c#3176).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#3176).z)), int32(0))>int32(0))) or (max(int32((int3d($c#3176).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0)))) then
      $Fluid#3173[$c#3176].velocityGradientX = vs_div_double_3(vs_mul_double_3(vv_sub_double_3($Fluid#3173[(($c#3176+{1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity, $Fluid#3173[(($c#3176+{-1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity), double(0.5)), double(0.0015841600529101))
      $Fluid#3173[$c#3176].velocityGradientY = vs_div_double_3(vs_mul_double_3(vv_sub_double_3($Fluid#3173[(($c#3176+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity, $Fluid#3173[(($c#3176+{0, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity), double(0.5)), double(0.0010188725631769))
      $Fluid#3173[$c#3176].velocityGradientZ = vs_div_double_3(vs_mul_double_3(vv_sub_double_3($Fluid#3173[(($c#3176+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity, $Fluid#3173[(($c#3176+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity), double(0.5)), double(0.0011352949453552))
    else
    end
  end
end
task Flow_UpdateAuxiliaryThermodynamics($dom#3302 : region#578(ispace#578(int3d), Fluid_columns), $Fluid#3304 : region#579(ispace#579(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#3304.pressure), writes($Fluid#3304.pressure), reads($Fluid#3304.rho), reads($Fluid#3304.rhoEnergy), reads($Fluid#3304.temperature), writes($Fluid#3304.temperature), reads($Fluid#3304.velocity), $dom#3302 <= $Fluid#3304
do
  for $c#3319 : int3d(Fluid_columns, $dom#3302) in $dom#3302 do
    if (not ((((((max(int32((uint64(int32(0))-int3d($c#3319).x)), int32(0))>int32(0)) or (max(int32((int3d($c#3319).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#3319).y)), int32(0))>int32(0))) or (max(int32((int3d($c#3319).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#3319).z)), int32(0))>int32(0))) or (max(int32((int3d($c#3319).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0)))) then
      var $kineticEnergy#3320 : double = ((double(0.5)*$Fluid#3304[$c#3319].rho)*dot_double_3($Fluid#3304[$c#3319].velocity, $Fluid#3304[$c#3319].velocity))
      var $pressure#3321 : double = ((double(0.8198141)-double(1))*($Fluid#3304[$c#3319].rhoEnergy-$kineticEnergy#3320))
      $Fluid#3304[$c#3319].pressure = $pressure#3321
      $Fluid#3304[$c#3319].temperature = ($pressure#3321/(double(0.7272462)*$Fluid#3304[$c#3319].rho))
    else
    end
  end
end
task Flow_UpdateGhostThermodynamicsStep1($dom#3355 : region#593(ispace#593(int3d), Fluid_columns), $Fluid#3357 : region#594(ispace#594(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#3357.pressure), reads($Fluid#3357.pressureBoundary), writes($Fluid#3357.pressureBoundary), reads($Fluid#3357.temperature), reads($Fluid#3357.temperatureBoundary), writes($Fluid#3357.temperatureBoundary), $dom#3355 <= $Fluid#3357
do
  for $c#3570 : int3d(Fluid_columns, $dom#3355) in $dom#3355 do
    if (max(int32((uint64(int32(0))-int3d($c#3570).x)), int32(0))>int32(0)) then
      var $c_bnd#3571 : int3d(Fluid_columns, $dom#3355) = $c#3570
      var $c_int#3572 : int3d = (($c#3570+{1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $bnd_temperature#3573 : double = double(-1)
      var $temp_wall#3574 : double = double(double(0))
      var $temperature#3575 : double = double(double(0))
      $temp_wall#3574 = $Fluid#3357[$c_int#3572].temperature
      if ($bnd_temperature#3573>double(0)) then
        $temp_wall#3574 = $bnd_temperature#3573
      else
      end
      $temperature#3575 = ((double(2)*$temp_wall#3574)-$Fluid#3357[$c_int#3572].temperature)
      $Fluid#3357[$c_bnd#3571].pressureBoundary = $Fluid#3357[$c_int#3572].pressure
      $Fluid#3357[$c_bnd#3571].temperatureBoundary = $temperature#3575
    else
    end
    if (max(int32((int3d($c#3570).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0)) then
      var $c_bnd#3576 : int3d(Fluid_columns, $dom#3355) = $c#3570
      var $c_int#3577 : int3d = (($c#3570+{-1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $bnd_temperature#3578 : double = double(-1)
      var $temp_wall#3579 : double = double(double(0))
      var $temperature#3580 : double = double(double(0))
      $temp_wall#3579 = $Fluid#3357[$c_int#3577].temperature
      if ($bnd_temperature#3578>double(0)) then
        $temp_wall#3579 = $bnd_temperature#3578
      else
      end
      $temperature#3580 = ((double(2)*$temp_wall#3579)-$Fluid#3357[$c_int#3577].temperature)
      $Fluid#3357[$c_bnd#3576].pressureBoundary = $Fluid#3357[$c_int#3577].pressure
      $Fluid#3357[$c_bnd#3576].temperatureBoundary = $temperature#3580
    else
    end
    if (max(int32((uint64(int32(1))-int3d($c#3570).y)), int32(0))>int32(0)) then
      var $c_bnd#3581 : int3d(Fluid_columns, $dom#3355) = $c#3570
      var $c_int#3582 : int3d = (($c#3570+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $bnd_temperature#3583 : double = double(-1)
      var $temp_wall#3584 : double = double(double(0))
      var $temperature#3585 : double = double(double(0))
      $temp_wall#3584 = $Fluid#3357[$c_int#3582].temperature
      if ($bnd_temperature#3583>double(0)) then
        $temp_wall#3584 = $bnd_temperature#3583
      else
      end
      $temperature#3585 = ((double(2)*$temp_wall#3584)-$Fluid#3357[$c_int#3582].temperature)
      $Fluid#3357[$c_bnd#3581].pressureBoundary = $Fluid#3357[$c_int#3582].pressure
      $Fluid#3357[$c_bnd#3581].temperatureBoundary = $temperature#3585
    else
    end
    if (max(int32((int3d($c#3570).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0)) then
      var $c_bnd#3586 : int3d(Fluid_columns, $dom#3355) = $c#3570
      var $c_int#3587 : int3d = (($c#3570+{0, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $bnd_temperature#3588 : double = double(-1)
      var $temp_wall#3589 : double = double(double(0))
      var $temperature#3590 : double = double(double(0))
      $temp_wall#3589 = $Fluid#3357[$c_int#3587].temperature
      if ($bnd_temperature#3588>double(0)) then
        $temp_wall#3589 = $bnd_temperature#3588
      else
      end
      $temperature#3590 = ((double(2)*$temp_wall#3589)-$Fluid#3357[$c_int#3587].temperature)
      $Fluid#3357[$c_bnd#3586].pressureBoundary = $Fluid#3357[$c_int#3587].pressure
      $Fluid#3357[$c_bnd#3586].temperatureBoundary = $temperature#3590
    else
    end
    if (max(int32((uint64(int32(1))-int3d($c#3570).z)), int32(0))>int32(0)) then
      var $c_bnd#3591 : int3d(Fluid_columns, $dom#3355) = $c#3570
      var $c_int#3592 : int3d = (($c#3570+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $bnd_temperature#3593 : double = double(-1)
      var $temp_wall#3594 : double = double(double(0))
      var $temperature#3595 : double = double(double(0))
      $temp_wall#3594 = $Fluid#3357[$c_int#3592].temperature
      if ($bnd_temperature#3593>double(0)) then
        $temp_wall#3594 = $bnd_temperature#3593
      else
      end
      $temperature#3595 = ((double(2)*$temp_wall#3594)-$Fluid#3357[$c_int#3592].temperature)
      $Fluid#3357[$c_bnd#3591].pressureBoundary = $Fluid#3357[$c_int#3592].pressure
      $Fluid#3357[$c_bnd#3591].temperatureBoundary = $temperature#3595
    else
    end
    if (max(int32((int3d($c#3570).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0)) then
      var $c_bnd#3596 : int3d(Fluid_columns, $dom#3355) = $c#3570
      var $c_int#3597 : int3d = (($c#3570+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $bnd_temperature#3598 : double = double(-1)
      var $temp_wall#3599 : double = double(double(0))
      var $temperature#3600 : double = double(double(0))
      $temp_wall#3599 = $Fluid#3357[$c_int#3597].temperature
      if ($bnd_temperature#3598>double(0)) then
        $temp_wall#3599 = $bnd_temperature#3598
      else
      end
      $temperature#3600 = ((double(2)*$temp_wall#3599)-$Fluid#3357[$c_int#3597].temperature)
      $Fluid#3357[$c_bnd#3596].pressureBoundary = $Fluid#3357[$c_int#3597].pressure
      $Fluid#3357[$c_bnd#3596].temperatureBoundary = $temperature#3600
    else
    end
  end
end
task Flow_UpdateGhostThermodynamicsStep2($dom#3816 : region#712(ispace#712(int3d), Fluid_columns), $Fluid#3818 : region#713(ispace#713(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#3818.pressure), writes($Fluid#3818.pressure), reads($Fluid#3818.pressureBoundary), reads($Fluid#3818.temperature), writes($Fluid#3818.temperature), reads($Fluid#3818.temperatureBoundary), $dom#3816 <= $Fluid#3818
do
  for $c#3821 : int3d(Fluid_columns, $dom#3816) in $dom#3816 do
    if ((((((max(int32((uint64(int32(0))-int3d($c#3821).x)), int32(0))>int32(0)) or (max(int32((int3d($c#3821).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#3821).y)), int32(0))>int32(0))) or (max(int32((int3d($c#3821).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#3821).z)), int32(0))>int32(0))) or (max(int32((int3d($c#3821).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0))) then
      $Fluid#3818[$c#3821].pressure = $Fluid#3818[$c#3821].pressureBoundary
      $Fluid#3818[$c#3821].temperature = $Fluid#3818[$c#3821].temperatureBoundary
    else
    end
  end
end
task Flow_UpdateGhostFieldsStep1($dom#3845 : region#723(ispace#723(int3d), Fluid_columns), $Fluid#3847 : region#724(ispace#724(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#3847.pressure), reads($Fluid#3847.pressureBoundary), writes($Fluid#3847.pressureBoundary), reads($Fluid#3847.rho), reads($Fluid#3847.rhoBoundary), writes($Fluid#3847.rhoBoundary), reads($Fluid#3847.rhoEnergyBoundary), writes($Fluid#3847.rhoEnergyBoundary), reads($Fluid#3847.rhoVelocity), reads($Fluid#3847.rhoVelocityBoundary), writes($Fluid#3847.rhoVelocityBoundary), reads($Fluid#3847.temperature), reads($Fluid#3847.temperatureBoundary), writes($Fluid#3847.temperatureBoundary), reads($Fluid#3847.velocityBoundary), writes($Fluid#3847.velocityBoundary), $dom#3845 <= $Fluid#3847
do
  for $c#4270 : int3d(Fluid_columns, $dom#3845) in $dom#3845 do
    if (max(int32((uint64(int32(0))-int3d($c#4270).x)), int32(0))>int32(0)) then
      var $c_bnd#4271 : int3d(Fluid_columns, $dom#3845) = $c#4270
      var $c_int#4272 : int3d = (($c#4270+{1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#4273 : double[3] = array(double(1), double(1), double(1))
      var $bnd_velocity#4274 : double[3] = array(double(0), double(0), double(0))
      var $bnd_temperature#4275 : double = double(-1)
      var $rho#4276 : double = double(double(0))
      var $temp_wall#4277 : double = double(double(0))
      var $temperature#4278 : double = double(double(0))
      var $velocity#4279 : double[3] = double[3](array(double(0), double(0), double(0)))
      var $cv#4280 : double = (double(0.7272462)/(double(0.8198141)-double(1)))
      $velocity#4279 = vv_add_double_3(vv_mul_double_3(vs_div_double_3($Fluid#3847[$c_int#4272].rhoVelocity, $Fluid#3847[$c_int#4272].rho), $sign#4273), $bnd_velocity#4274)
      $temp_wall#4277 = $Fluid#3847[$c_int#4272].temperature
      if ($bnd_temperature#4275>double(0)) then
        $temp_wall#4277 = $bnd_temperature#4275
      else
      end
      $temperature#4278 = ((double(2)*$temp_wall#4277)-$Fluid#3847[$c_int#4272].temperature)
      $rho#4276 = ($Fluid#3847[$c_int#4272].pressure/(double(0.7272462)*$temperature#4278))
      $Fluid#3847[$c_bnd#4271].rhoBoundary = $rho#4276
      $Fluid#3847[$c_bnd#4271].rhoVelocityBoundary = vs_mul_double_3($velocity#4279, $rho#4276)
      $Fluid#3847[$c_bnd#4271].rhoEnergyBoundary = ($rho#4276*(($cv#4280*$temperature#4278)+(double(0.5)*dot_double_3($velocity#4279, $velocity#4279))))
      $Fluid#3847[$c_bnd#4271].velocityBoundary = $velocity#4279
      $Fluid#3847[$c_bnd#4271].pressureBoundary = $Fluid#3847[$c_int#4272].pressure
      $Fluid#3847[$c_bnd#4271].temperatureBoundary = $temperature#4278
    else
    end
    if (max(int32((int3d($c#4270).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0)) then
      var $c_bnd#4281 : int3d(Fluid_columns, $dom#3845) = $c#4270
      var $c_int#4282 : int3d = (($c#4270+{-1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#4283 : double[3] = array(double(1), double(1), double(1))
      var $bnd_velocity#4284 : double[3] = array(double(0), double(0), double(0))
      var $bnd_temperature#4285 : double = double(-1)
      var $rho#4286 : double = double(double(0))
      var $temp_wall#4287 : double = double(double(0))
      var $temperature#4288 : double = double(double(0))
      var $velocity#4289 : double[3] = double[3](array(double(0), double(0), double(0)))
      var $cv#4290 : double = (double(0.7272462)/(double(0.8198141)-double(1)))
      $velocity#4289 = vv_add_double_3(vv_mul_double_3(vs_div_double_3($Fluid#3847[$c_int#4282].rhoVelocity, $Fluid#3847[$c_int#4282].rho), $sign#4283), $bnd_velocity#4284)
      $temp_wall#4287 = $Fluid#3847[$c_int#4282].temperature
      if ($bnd_temperature#4285>double(0)) then
        $temp_wall#4287 = $bnd_temperature#4285
      else
      end
      $temperature#4288 = ((double(2)*$temp_wall#4287)-$Fluid#3847[$c_int#4282].temperature)
      $rho#4286 = ($Fluid#3847[$c_int#4282].pressure/(double(0.7272462)*$temperature#4288))
      $Fluid#3847[$c_bnd#4281].rhoBoundary = $rho#4286
      $Fluid#3847[$c_bnd#4281].rhoVelocityBoundary = vs_mul_double_3($velocity#4289, $rho#4286)
      $Fluid#3847[$c_bnd#4281].rhoEnergyBoundary = ($rho#4286*(($cv#4290*$temperature#4288)+(double(0.5)*dot_double_3($velocity#4289, $velocity#4289))))
      $Fluid#3847[$c_bnd#4281].velocityBoundary = $velocity#4289
      $Fluid#3847[$c_bnd#4281].pressureBoundary = $Fluid#3847[$c_int#4282].pressure
      $Fluid#3847[$c_bnd#4281].temperatureBoundary = $temperature#4288
    else
    end
    if (max(int32((uint64(int32(1))-int3d($c#4270).y)), int32(0))>int32(0)) then
      var $c_bnd#4291 : int3d(Fluid_columns, $dom#3845) = $c#4270
      var $c_int#4292 : int3d = (($c#4270+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#4293 : double[3] = array(double(-1), double(-1), double(-1))
      var $bnd_velocity#4294 : double[3] = array(double(1.9836594), double(1.7248182), double(0.6582188))
      var $bnd_temperature#4295 : double = double(-1)
      var $rho#4296 : double = double(double(0))
      var $temp_wall#4297 : double = double(double(0))
      var $temperature#4298 : double = double(double(0))
      var $velocity#4299 : double[3] = double[3](array(double(0), double(0), double(0)))
      var $cv#4300 : double = (double(0.7272462)/(double(0.8198141)-double(1)))
      $velocity#4299 = vv_add_double_3(vv_mul_double_3(vs_div_double_3($Fluid#3847[$c_int#4292].rhoVelocity, $Fluid#3847[$c_int#4292].rho), $sign#4293), $bnd_velocity#4294)
      $temp_wall#4297 = $Fluid#3847[$c_int#4292].temperature
      if ($bnd_temperature#4295>double(0)) then
        $temp_wall#4297 = $bnd_temperature#4295
      else
      end
      $temperature#4298 = ((double(2)*$temp_wall#4297)-$Fluid#3847[$c_int#4292].temperature)
      $rho#4296 = ($Fluid#3847[$c_int#4292].pressure/(double(0.7272462)*$temperature#4298))
      $Fluid#3847[$c_bnd#4291].rhoBoundary = $rho#4296
      $Fluid#3847[$c_bnd#4291].rhoVelocityBoundary = vs_mul_double_3($velocity#4299, $rho#4296)
      $Fluid#3847[$c_bnd#4291].rhoEnergyBoundary = ($rho#4296*(($cv#4300*$temperature#4298)+(double(0.5)*dot_double_3($velocity#4299, $velocity#4299))))
      $Fluid#3847[$c_bnd#4291].velocityBoundary = $velocity#4299
      $Fluid#3847[$c_bnd#4291].pressureBoundary = $Fluid#3847[$c_int#4292].pressure
      $Fluid#3847[$c_bnd#4291].temperatureBoundary = $temperature#4298
    else
    end
    if (max(int32((int3d($c#4270).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0)) then
      var $c_bnd#4301 : int3d(Fluid_columns, $dom#3845) = $c#4270
      var $c_int#4302 : int3d = (($c#4270+{0, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#4303 : double[3] = array(double(-1), double(-1), double(-1))
      var $bnd_velocity#4304 : double[3] = array(double(0.8413202), double(1.3283396), double(1.134841))
      var $bnd_temperature#4305 : double = double(-1)
      var $rho#4306 : double = double(double(0))
      var $temp_wall#4307 : double = double(double(0))
      var $temperature#4308 : double = double(double(0))
      var $velocity#4309 : double[3] = double[3](array(double(0), double(0), double(0)))
      var $cv#4310 : double = (double(0.7272462)/(double(0.8198141)-double(1)))
      $velocity#4309 = vv_add_double_3(vv_mul_double_3(vs_div_double_3($Fluid#3847[$c_int#4302].rhoVelocity, $Fluid#3847[$c_int#4302].rho), $sign#4303), $bnd_velocity#4304)
      $temp_wall#4307 = $Fluid#3847[$c_int#4302].temperature
      if ($bnd_temperature#4305>double(0)) then
        $temp_wall#4307 = $bnd_temperature#4305
      else
      end
      $temperature#4308 = ((double(2)*$temp_wall#4307)-$Fluid#3847[$c_int#4302].temperature)
      $rho#4306 = ($Fluid#3847[$c_int#4302].pressure/(double(0.7272462)*$temperature#4308))
      $Fluid#3847[$c_bnd#4301].rhoBoundary = $rho#4306
      $Fluid#3847[$c_bnd#4301].rhoVelocityBoundary = vs_mul_double_3($velocity#4309, $rho#4306)
      $Fluid#3847[$c_bnd#4301].rhoEnergyBoundary = ($rho#4306*(($cv#4310*$temperature#4308)+(double(0.5)*dot_double_3($velocity#4309, $velocity#4309))))
      $Fluid#3847[$c_bnd#4301].velocityBoundary = $velocity#4309
      $Fluid#3847[$c_bnd#4301].pressureBoundary = $Fluid#3847[$c_int#4302].pressure
      $Fluid#3847[$c_bnd#4301].temperatureBoundary = $temperature#4308
    else
    end
    if (max(int32((uint64(int32(1))-int3d($c#4270).z)), int32(0))>int32(0)) then
      var $c_bnd#4311 : int3d(Fluid_columns, $dom#3845) = $c#4270
      var $c_int#4312 : int3d = (($c#4270+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#4313 : double[3] = array(double(1), double(1), double(-1))
      var $bnd_velocity#4314 : double[3] = array(double(0), double(0), double(0))
      var $bnd_temperature#4315 : double = double(-1)
      var $rho#4316 : double = double(double(0))
      var $temp_wall#4317 : double = double(double(0))
      var $temperature#4318 : double = double(double(0))
      var $velocity#4319 : double[3] = double[3](array(double(0), double(0), double(0)))
      var $cv#4320 : double = (double(0.7272462)/(double(0.8198141)-double(1)))
      $velocity#4319 = vv_add_double_3(vv_mul_double_3(vs_div_double_3($Fluid#3847[$c_int#4312].rhoVelocity, $Fluid#3847[$c_int#4312].rho), $sign#4313), $bnd_velocity#4314)
      $temp_wall#4317 = $Fluid#3847[$c_int#4312].temperature
      if ($bnd_temperature#4315>double(0)) then
        $temp_wall#4317 = $bnd_temperature#4315
      else
      end
      $temperature#4318 = ((double(2)*$temp_wall#4317)-$Fluid#3847[$c_int#4312].temperature)
      $rho#4316 = ($Fluid#3847[$c_int#4312].pressure/(double(0.7272462)*$temperature#4318))
      $Fluid#3847[$c_bnd#4311].rhoBoundary = $rho#4316
      $Fluid#3847[$c_bnd#4311].rhoVelocityBoundary = vs_mul_double_3($velocity#4319, $rho#4316)
      $Fluid#3847[$c_bnd#4311].rhoEnergyBoundary = ($rho#4316*(($cv#4320*$temperature#4318)+(double(0.5)*dot_double_3($velocity#4319, $velocity#4319))))
      $Fluid#3847[$c_bnd#4311].velocityBoundary = $velocity#4319
      $Fluid#3847[$c_bnd#4311].pressureBoundary = $Fluid#3847[$c_int#4312].pressure
      $Fluid#3847[$c_bnd#4311].temperatureBoundary = $temperature#4318
    else
    end
    if (max(int32((int3d($c#4270).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0)) then
      var $c_bnd#4321 : int3d(Fluid_columns, $dom#3845) = $c#4270
      var $c_int#4322 : int3d = (($c#4270+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#4323 : double[3] = array(double(1), double(1), double(-1))
      var $bnd_velocity#4324 : double[3] = array(double(0), double(0), double(0))
      var $bnd_temperature#4325 : double = double(-1)
      var $rho#4326 : double = double(double(0))
      var $temp_wall#4327 : double = double(double(0))
      var $temperature#4328 : double = double(double(0))
      var $velocity#4329 : double[3] = double[3](array(double(0), double(0), double(0)))
      var $cv#4330 : double = (double(0.7272462)/(double(0.8198141)-double(1)))
      $velocity#4329 = vv_add_double_3(vv_mul_double_3(vs_div_double_3($Fluid#3847[$c_int#4322].rhoVelocity, $Fluid#3847[$c_int#4322].rho), $sign#4323), $bnd_velocity#4324)
      $temp_wall#4327 = $Fluid#3847[$c_int#4322].temperature
      if ($bnd_temperature#4325>double(0)) then
        $temp_wall#4327 = $bnd_temperature#4325
      else
      end
      $temperature#4328 = ((double(2)*$temp_wall#4327)-$Fluid#3847[$c_int#4322].temperature)
      $rho#4326 = ($Fluid#3847[$c_int#4322].pressure/(double(0.7272462)*$temperature#4328))
      $Fluid#3847[$c_bnd#4321].rhoBoundary = $rho#4326
      $Fluid#3847[$c_bnd#4321].rhoVelocityBoundary = vs_mul_double_3($velocity#4329, $rho#4326)
      $Fluid#3847[$c_bnd#4321].rhoEnergyBoundary = ($rho#4326*(($cv#4330*$temperature#4328)+(double(0.5)*dot_double_3($velocity#4329, $velocity#4329))))
      $Fluid#3847[$c_bnd#4321].velocityBoundary = $velocity#4329
      $Fluid#3847[$c_bnd#4321].pressureBoundary = $Fluid#3847[$c_int#4322].pressure
      $Fluid#3847[$c_bnd#4321].temperatureBoundary = $temperature#4328
    else
    end
  end
end
task Flow_UpdateGhostFieldsStep2($dom#4900 : region#938(ispace#938(int3d), Fluid_columns), $Fluid#4902 : region#939(ispace#939(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#4902.pressure), writes($Fluid#4902.pressure), reads($Fluid#4902.pressureBoundary), reads($Fluid#4902.rho), writes($Fluid#4902.rho), reads($Fluid#4902.rhoBoundary), reads($Fluid#4902.rhoEnergy), writes($Fluid#4902.rhoEnergy), reads($Fluid#4902.rhoEnergyBoundary), reads($Fluid#4902.rhoVelocity), writes($Fluid#4902.rhoVelocity), reads($Fluid#4902.rhoVelocityBoundary), reads($Fluid#4902.temperature), writes($Fluid#4902.temperature), reads($Fluid#4902.temperatureBoundary), $dom#4900 <= $Fluid#4902
do
  for $c#4905 : int3d(Fluid_columns, $dom#4900) in $dom#4900 do
    if ((((((max(int32((uint64(int32(0))-int3d($c#4905).x)), int32(0))>int32(0)) or (max(int32((int3d($c#4905).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#4905).y)), int32(0))>int32(0))) or (max(int32((int3d($c#4905).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#4905).z)), int32(0))>int32(0))) or (max(int32((int3d($c#4905).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0))) then
      $Fluid#4902[$c#4905].rho = $Fluid#4902[$c#4905].rhoBoundary
      $Fluid#4902[$c#4905].rhoVelocity = $Fluid#4902[$c#4905].rhoVelocityBoundary
      $Fluid#4902[$c#4905].rhoEnergy = $Fluid#4902[$c#4905].rhoEnergyBoundary
      $Fluid#4902[$c#4905].pressure = $Fluid#4902[$c#4905].pressureBoundary
      $Fluid#4902[$c#4905].temperature = $Fluid#4902[$c#4905].temperatureBoundary
    else
    end
  end
end
task numberOfInteriorCells($dom#4929 : region#949(ispace#949(int3d), Fluid_columns), $Fluid#4931 : region#950(ispace#950(int3d), Fluid_columns)) : int64
-- leaf (false), inner (false), idempotent (false)
where
  $dom#4929 <= $Fluid#4931
do
  var $acc#4936 : int64 = int64(0)
  for $c#4937 : int3d(Fluid_columns, $dom#4929) in $dom#4929 do
    if (not ((((((max(int32((uint64(int32(0))-int3d($c#4937).x)), int32(0))>int32(0)) or (max(int32((int3d($c#4937).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#4937).y)), int32(0))>int32(0))) or (max(int32((int3d($c#4937).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#4937).z)), int32(0))>int32(0))) or (max(int32((int3d($c#4937).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0)))) then
      $acc#4936 += int64(int32(1))
    else
    end
  end
  return $acc#4936
end
task areaInterior($dom#4964 : region#962(ispace#962(int3d), Fluid_columns), $Fluid#4966 : region#963(ispace#963(int3d), Fluid_columns)) : double
-- leaf (false), inner (false), idempotent (false)
where
  $dom#4964 <= $Fluid#4966
do
  var $acc#4971 : double = double(0)
  for $c#4972 : int3d(Fluid_columns, $dom#4964) in $dom#4964 do
    if (not ((((((max(int32((uint64(int32(0))-int3d($c#4972).x)), int32(0))>int32(0)) or (max(int32((int3d($c#4972).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#4972).y)), int32(0))>int32(0))) or (max(int32((int3d($c#4972).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#4972).z)), int32(0))>int32(0))) or (max(int32((int3d($c#4972).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0)))) then
      $acc#4971 += double(1.8324309961038e-09)
    else
    end
  end
  return $acc#4971
end
task averagePressure($dom#4999 : region#975(ispace#975(int3d), Fluid_columns), $Fluid#5001 : region#976(ispace#976(int3d), Fluid_columns)) : double
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#5001.pressure), $dom#4999 <= $Fluid#5001
do
  var $acc#5006 : double = double(0)
  for $c#5007 : int3d(Fluid_columns, $dom#4999) in $dom#4999 do
    if (not ((((((max(int32((uint64(int32(0))-int3d($c#5007).x)), int32(0))>int32(0)) or (max(int32((int3d($c#5007).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#5007).y)), int32(0))>int32(0))) or (max(int32((int3d($c#5007).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#5007).z)), int32(0))>int32(0))) or (max(int32((int3d($c#5007).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0)))) then
      $acc#5006 += ($Fluid#5001[$c#5007].pressure*double(1.8324309961038e-09))
    else
    end
  end
  return $acc#5006
end
task averageTemperature($dom#5034 : region#988(ispace#988(int3d), Fluid_columns), $Fluid#5036 : region#989(ispace#989(int3d), Fluid_columns)) : double
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#5036.temperature), $dom#5034 <= $Fluid#5036
do
  var $acc#5041 : double = double(0)
  for $c#5042 : int3d(Fluid_columns, $dom#5034) in $dom#5034 do
    if (not ((((((max(int32((uint64(int32(0))-int3d($c#5042).x)), int32(0))>int32(0)) or (max(int32((int3d($c#5042).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#5042).y)), int32(0))>int32(0))) or (max(int32((int3d($c#5042).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#5042).z)), int32(0))>int32(0))) or (max(int32((int3d($c#5042).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0)))) then
      $acc#5041 += ($Fluid#5036[$c#5042].temperature*double(1.8324309961038e-09))
    else
    end
  end
  return $acc#5041
end
task averageKineticEnergy($dom#5069 : region#1001(ispace#1001(int3d), Fluid_columns), $Fluid#5071 : region#1002(ispace#1002(int3d), Fluid_columns)) : double
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#5071.kineticEnergy), $dom#5069 <= $Fluid#5071
do
  var $acc#5076 : double = double(0)
  for $c#5077 : int3d(Fluid_columns, $dom#5069) in $dom#5069 do
    if (not ((((((max(int32((uint64(int32(0))-int3d($c#5077).x)), int32(0))>int32(0)) or (max(int32((int3d($c#5077).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#5077).y)), int32(0))>int32(0))) or (max(int32((int3d($c#5077).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#5077).z)), int32(0))>int32(0))) or (max(int32((int3d($c#5077).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0)))) then
      $acc#5076 += ($Fluid#5071[$c#5077].kineticEnergy*double(1.8324309961038e-09))
    else
    end
  end
  return $acc#5076
end
task minTemperature($dom#5104 : region#1014(ispace#1014(int3d), Fluid_columns), $Fluid#5106 : region#1015(ispace#1015(int3d), Fluid_columns)) : double
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#5106.temperature), $dom#5104 <= $Fluid#5106
do
  var $acc#5111 : double = inf
  for $c#5112 : int3d(Fluid_columns, $dom#5104) in $dom#5104 do
    if (not ((((((max(int32((uint64(int32(0))-int3d($c#5112).x)), int32(0))>int32(0)) or (max(int32((int3d($c#5112).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#5112).y)), int32(0))>int32(0))) or (max(int32((int3d($c#5112).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#5112).z)), int32(0))>int32(0))) or (max(int32((int3d($c#5112).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0)))) then
      $acc#5111 min= $Fluid#5106[$c#5112].temperature
    else
    end
  end
  return $acc#5111
end
task maxTemperature($dom#5139 : region#1027(ispace#1027(int3d), Fluid_columns), $Fluid#5141 : region#1028(ispace#1028(int3d), Fluid_columns)) : double
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#5141.temperature), $dom#5139 <= $Fluid#5141
do
  var $acc#5146 : double = -inf
  for $c#5147 : int3d(Fluid_columns, $dom#5139) in $dom#5139 do
    if (not ((((((max(int32((uint64(int32(0))-int3d($c#5147).x)), int32(0))>int32(0)) or (max(int32((int3d($c#5147).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#5147).y)), int32(0))>int32(0))) or (max(int32((int3d($c#5147).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#5147).z)), int32(0))>int32(0))) or (max(int32((int3d($c#5147).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0)))) then
      $acc#5146 max= $Fluid#5141[$c#5147].temperature
    else
    end
  end
  return $acc#5146
end
task Particles_IntegrateQuantities($dom#5174 : region#1040(ispace#1040(int1d), particles_columns), $particles#5176 : region#1041(ispace#1041(int1d), particles_columns)) : double
-- leaf (false), inner (false), idempotent (false)
where
  reads($particles#5176.particle_temperature), reads($particles#5176.__valid), $dom#5174 <= $particles#5176
do
  var $acc#5181 : double = double(0)
  for $p#5182 : int1d(particles_columns, $dom#5174) in $dom#5174 do
    if $particles#5176[$p#5182].__valid then
      $acc#5181 += $particles#5176[$p#5182].particle_temperature
    else
    end
  end
  return $acc#5181
end
task print($5210 : double)
-- leaf (false), inner (false), idempotent (false)
  printf("\n Current time step: %2.6e s.\n", $5210)
end
task print_($5216 : double, $5217 : double)
-- leaf (false), inner (false), idempotent (false)
  printf(" Min Flow Temp: %11.6f K. Max Flow Temp: %11.6f K.\n", $5216, $5217)
end
task print__($5224 : int64)
-- leaf (false), inner (false), idempotent (false)
  printf(" Current number of particles: %d.\n", $5224)
end
task print___()
-- leaf (false), inner (false), idempotent (false)
  printf("\n")
end
task print____()
-- leaf (false), inner (false), idempotent (false)
  printf("    Iter     Time(s)   Avg Press    Avg Temp      Avg KE  Particle T\n")
end
task print_____($5239 : int32, $5240 : double, $5241 : double, $5242 : double, $5243 : double, $5244 : double)
-- leaf (false), inner (false), idempotent (false)
  printf("%8d %11.6f %11.6f %11.6f %11.6f %11.6f\n", $5239, $5240, $5241, $5242, $5243, $5244)
end
terra Fluid_hdf5create_rho_pressure_velocity(filename : &int8) : {}
    var fid : int32 = H5Fcreate(filename, [uint32](2), 0, 0)
    var dataSpace : int32
    var sizes : uint64[3]
    [&uint64](sizes)[2] = [uint64](378)
    [&uint64](sizes)[1] = [uint64](833)
    [&uint64](sizes)[0] = [uint64](734)
    dataSpace = H5Screate_simple(3, [&uint64](sizes), [&uint64](0))
    var hType : int32 = extern global H5T_IEEE_F64LE_g : int32
    var $dataSet : int32 = H5Dcreate2(fid, "rho", hType, dataSpace, 0, 0, 0)
    var hType$1 : int32 = extern global H5T_IEEE_F64LE_g : int32
    var $dataSet$1 : int32 = H5Dcreate2(fid, "pressure", hType$1, dataSpace, 0, 0, 0)
    var dims : uint64[1]
    [&uint64](dims)[0] = [uint64](3)
    var elemType : int32 = extern global H5T_IEEE_F64LE_g : int32
    var $arrayType : int32 = H5Tarray_create2(elemType, [uint32](1), [&uint64](dims))
    var hType$2 : int32 = $arrayType
    var $dataSet$2 : int32 = H5Dcreate2(fid, "velocity", hType$2, dataSpace, 0, 0, 0)
    H5Dclose($dataSet$2)
    H5Tclose($arrayType)
    H5Dclose($dataSet$1)
    H5Dclose($dataSet)
    H5Sclose(dataSpace)
    H5Fclose(fid)
end
terra particles_hdf5create_cell_position_particle_velocity_particle_temperature_diameter___valid(filename : &int8) : {}
    var fid : int32 = H5Fcreate(filename, [uint32](2), 0, 0)
    var dataSpace : int32
    var sizes : uint64[1]
    [&uint64](sizes)[0] = [uint64](15 * 24)
    dataSpace = H5Screate_simple(1, [&uint64](sizes), [&uint64](0))
    var $int3dType : int32 = H5Tcreate([uint32](6), [uint64](24))
    H5Tinsert($int3dType, "x", [uint64](0), extern global H5T_STD_I64LE_g : int32)
    H5Tinsert($int3dType, "y", [uint64](8), extern global H5T_STD_I64LE_g : int32)
    H5Tinsert($int3dType, "z", [uint64](16), extern global H5T_STD_I64LE_g : int32)
    var $dataSet : int32 = H5Dcreate2(fid, "cell", $int3dType, dataSpace, 0, 0, 0)
    var dims : uint64[1]
    [&uint64](dims)[0] = [uint64](3)
    var elemType : int32 = extern global H5T_IEEE_F64LE_g : int32
    var $arrayType : int32 = H5Tarray_create2(elemType, [uint32](1), [&uint64](dims))
    var hType : int32 = $arrayType
    var $dataSet$1 : int32 = H5Dcreate2(fid, "position", hType, dataSpace, 0, 0, 0)
    var dims$1 : uint64[1]
    [&uint64](dims$1)[0] = [uint64](3)
    var elemType$1 : int32 = extern global H5T_IEEE_F64LE_g : int32
    var $arrayType$1 : int32 = H5Tarray_create2(elemType$1, [uint32](1), [&uint64](dims$1))
    var hType$1 : int32 = $arrayType$1
    var $dataSet$2 : int32 = H5Dcreate2(fid, "particle_velocity", hType$1, dataSpace, 0, 0, 0)
    var hType$2 : int32 = extern global H5T_IEEE_F64LE_g : int32
    var $dataSet$3 : int32 = H5Dcreate2(fid, "particle_temperature", hType$2, dataSpace, 0, 0, 0)
    var hType$3 : int32 = extern global H5T_IEEE_F64LE_g : int32
    var $dataSet$4 : int32 = H5Dcreate2(fid, "diameter", hType$3, dataSpace, 0, 0, 0)
    var hType$4 : int32 = extern global H5T_STD_U8LE_g : int32
    var $dataSet$5 : int32 = H5Dcreate2(fid, "__valid", hType$4, dataSpace, 0, 0, 0)
    H5Dclose($dataSet$5)
    H5Dclose($dataSet$4)
    H5Dclose($dataSet$3)
    H5Dclose($dataSet$2)
    H5Tclose($arrayType$1)
    H5Dclose($dataSet$1)
    H5Tclose($arrayType)
    H5Dclose($dataSet)
    H5Tclose($int3dType)
    H5Sclose(dataSpace)
    H5Fclose(fid)
end
task GetSoundSpeed($temperature#5301 : double) : double
-- leaf (false), inner (false), idempotent (false)
  return [regentlib.sqrt(double)](((double(0.8198141)*double(0.7272462))*$temperature#5301))
end
task calculateConvectiveSpectralRadius($dom#5295 : region#1069(ispace#1069(int3d), Fluid_columns), $Fluid#5297 : region#1070(ispace#1070(int3d), Fluid_columns)) : double
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#5297.convectiveSpectralRadius), writes($Fluid#5297.convectiveSpectralRadius), reads($Fluid#5297.temperature), reads($Fluid#5297.velocity), $dom#5295 <= $Fluid#5297
do
  var $acc#5306 : double = -inf
  for $c#5307 : int3d(Fluid_columns, $dom#5295) in $dom#5295 do
    var $5320 : double
    do
      var $5319 : double = $Fluid#5297[$c#5307].temperature
      $5320 = [regentlib.sqrt(double)](((double(0.8198141)*double(0.7272462))*$5319))
    end
    var $dummy#5321 : int32 = 0
    $Fluid#5297[$c#5307].convectiveSpectralRadius = (((([regentlib.fabs(double)]($Fluid#5297[$c#5307].velocity[int32(0)])/double(0.0015841600529101))+([regentlib.fabs(double)]($Fluid#5297[$c#5307].velocity[int32(1)])/double(0.0010188725631769)))+([regentlib.fabs(double)]($Fluid#5297[$c#5307].velocity[int32(2)])/double(0.0011352949453552)))+($5320*[regentlib.sqrt(double)](double(2137631.4894183))))
    $acc#5306 max= $Fluid#5297[$c#5307].convectiveSpectralRadius
  end
  return $acc#5306
end
task GetDynamicViscosity($temperature#5387 : double) : double
-- leaf (false), inner (false), idempotent (false)
  var $viscosity#5391 : double = double(double(0))
  if (int32(3003)==int32(3001)) then
    $viscosity#5391 = double(0.2454194)
  else
    if (int32(3003)==int32(3002)) then
      $viscosity#5391 = (double(0.9683056)*pow(($temperature#5387/double(0.9913615)), double(0.75)))
    else
      if (int32(3003)==int32(3003)) then
        $viscosity#5391 = ((double(0.5815515)*pow(($temperature#5387/double(0.6935851)), (double(3)/double(2))))*((double(0.6935851)+double(0.2772065))/($temperature#5387+double(0.2772065))))
      else
        std.assert(false, "(Liszt assertion)")
      end
    end
  end
  return $viscosity#5391
end
task calculateViscousSpectralRadius($dom#5380 : region#1103(ispace#1103(int3d), Fluid_columns), $Fluid#5382 : region#1104(ispace#1104(int3d), Fluid_columns)) : double
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#5382.rho), reads($Fluid#5382.sgsEddyViscosity), reads($Fluid#5382.temperature), reads($Fluid#5382.viscousSpectralRadius), writes($Fluid#5382.viscousSpectralRadius), $dom#5380 <= $Fluid#5382
do
  var $acc#5415 : double = -inf
  for $c#5416 : int3d(Fluid_columns, $dom#5380) in $dom#5380 do
    var $5425 : double
    do
      var $5424 : double = $Fluid#5382[$c#5416].temperature
      var $viscosity#5391 : double = double(double(0))
      if (int32(3003)==int32(3001)) then
        $viscosity#5391 = double(0.2454194)
      else
        if (int32(3003)==int32(3002)) then
          $viscosity#5391 = (double(0.9683056)*pow(($5424/double(0.9913615)), double(0.75)))
        else
          if (int32(3003)==int32(3003)) then
            $viscosity#5391 = ((double(0.5815515)*pow(($5424/double(0.6935851)), (double(3)/double(2))))*((double(0.6935851)+double(0.2772065))/($5424+double(0.2772065))))
          else
            std.assert(false, "(Liszt assertion)")
          end
        end
      end
      $5425 = $viscosity#5391
    end
    var $dummy#5426 : int32 = 0
    var $dynamicViscosity#5417 : double = $5425
    var $eddyViscosity#5418 : double = $Fluid#5382[$c#5416].sgsEddyViscosity
    $Fluid#5382[$c#5416].viscousSpectralRadius = ((((double(2)*($dynamicViscosity#5417+$eddyViscosity#5418))/$Fluid#5382[$c#5416].rho)*double(2137631.4894183))*double(4))
    $acc#5415 max= $Fluid#5382[$c#5416].viscousSpectralRadius
  end
  return $acc#5415
end
task calculateHeatConductionSpectralRadius($dom#5471 : region#1134(ispace#1134(int3d), Fluid_columns), $Fluid#5473 : region#1135(ispace#1135(int3d), Fluid_columns)) : double
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#5473.heatConductionSpectralRadius), writes($Fluid#5473.heatConductionSpectralRadius), reads($Fluid#5473.rho), reads($Fluid#5473.sgsEddyKappa), reads($Fluid#5473.temperature), $dom#5471 <= $Fluid#5473
do
  var $acc#5494 : double = -inf
  for $c#5495 : int3d(Fluid_columns, $dom#5471) in $dom#5471 do
    var $5508 : double
    do
      var $5507 : double = $Fluid#5473[$c#5495].temperature
      var $viscosity#5391 : double = double(double(0))
      if (int32(3003)==int32(3001)) then
        $viscosity#5391 = double(0.2454194)
      else
        if (int32(3003)==int32(3002)) then
          $viscosity#5391 = (double(0.9683056)*pow(($5507/double(0.9913615)), double(0.75)))
        else
          if (int32(3003)==int32(3003)) then
            $viscosity#5391 = ((double(0.5815515)*pow(($5507/double(0.6935851)), (double(3)/double(2))))*((double(0.6935851)+double(0.2772065))/($5507+double(0.2772065))))
          else
            std.assert(false, "(Liszt assertion)")
          end
        end
      end
      $5508 = $viscosity#5391
    end
    var $dummy#5509 : int32 = 0
    var $dynamicViscosity#5496 : double = $5508
    var $cv#5497 : double = (double(0.7272462)/(double(0.8198141)-double(1)))
    var $cp#5498 : double = (double(0.8198141)*$cv#5497)
    var $kappa#5499 : double = (($cp#5498/double(0.1565651))*$dynamicViscosity#5496)
    $Fluid#5473[$c#5495].heatConductionSpectralRadius = (((($kappa#5499+$Fluid#5473[$c#5495].sgsEddyKappa)/($cv#5497*$Fluid#5473[$c#5495].rho))*double(2137631.4894183))*double(4))
    $acc#5494 max= $Fluid#5473[$c#5495].heatConductionSpectralRadius
  end
  return $acc#5494
end
task Flow_InitializeTemporaries($dom#5558 : region#1167(ispace#1167(int3d), Fluid_columns), $Fluid#5560 : region#1168(ispace#1168(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#5560.rho), reads($Fluid#5560.rhoEnergy), reads($Fluid#5560.rhoEnergy_new), writes($Fluid#5560.rhoEnergy_new), reads($Fluid#5560.rhoEnergy_old), writes($Fluid#5560.rhoEnergy_old), reads($Fluid#5560.rhoVelocity), reads($Fluid#5560.rhoVelocity_new), writes($Fluid#5560.rhoVelocity_new), reads($Fluid#5560.rhoVelocity_old), writes($Fluid#5560.rhoVelocity_old), reads($Fluid#5560.rho_new), writes($Fluid#5560.rho_new), reads($Fluid#5560.rho_old), writes($Fluid#5560.rho_old), $dom#5558 <= $Fluid#5560
do
  for $c#5563 : int3d(Fluid_columns, $dom#5558) in $dom#5558 do
    $Fluid#5560[$c#5563].rho_old = $Fluid#5560[$c#5563].rho
    $Fluid#5560[$c#5563].rhoVelocity_old = $Fluid#5560[$c#5563].rhoVelocity
    $Fluid#5560[$c#5563].rhoEnergy_old = $Fluid#5560[$c#5563].rhoEnergy
    $Fluid#5560[$c#5563].rho_new = $Fluid#5560[$c#5563].rho
    $Fluid#5560[$c#5563].rhoVelocity_new = $Fluid#5560[$c#5563].rhoVelocity
    $Fluid#5560[$c#5563].rhoEnergy_new = $Fluid#5560[$c#5563].rhoEnergy
  end
end
task Particles_InitializeTemporaries($dom#5587 : region#1178(ispace#1178(int1d), particles_columns), $particles#5589 : region#1179(ispace#1179(int1d), particles_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($particles#5589.particle_temperature), reads($particles#5589.particle_velocity), reads($particles#5589.position), reads($particles#5589.position_new), writes($particles#5589.position_new), reads($particles#5589.position_old), writes($particles#5589.position_old), reads($particles#5589.temperature_new), writes($particles#5589.temperature_new), reads($particles#5589.temperature_old), writes($particles#5589.temperature_old), reads($particles#5589.velocity_new), writes($particles#5589.velocity_new), reads($particles#5589.velocity_old), writes($particles#5589.velocity_old), reads($particles#5589.__valid), $dom#5587 <= $particles#5589
do
  for $p#5592 : int1d(particles_columns, $dom#5587) in $dom#5587 do
    if $particles#5589[$p#5592].__valid then
      $particles#5589[$p#5592].position_old = $particles#5589[$p#5592].position
      $particles#5589[$p#5592].velocity_old = $particles#5589[$p#5592].particle_velocity
      $particles#5589[$p#5592].temperature_old = $particles#5589[$p#5592].particle_temperature
      $particles#5589[$p#5592].position_new = $particles#5589[$p#5592].position
      $particles#5589[$p#5592].velocity_new = $particles#5589[$p#5592].particle_velocity
      $particles#5589[$p#5592].temperature_new = $particles#5589[$p#5592].particle_temperature
    else
    end
  end
end
task Flow_InitializeTimeDerivatives($dom#5616 : region#1189(ispace#1189(int3d), Fluid_columns), $Fluid#5618 : region#1190(ispace#1190(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#5618.pressure), reads($Fluid#5618.rhoEnergy), reads($Fluid#5618.rhoEnergy_t), writes($Fluid#5618.rhoEnergy_t), reads($Fluid#5618.rhoEnthalpy), writes($Fluid#5618.rhoEnthalpy), reads($Fluid#5618.rhoVelocity_t), writes($Fluid#5618.rhoVelocity_t), reads($Fluid#5618.rho_t), writes($Fluid#5618.rho_t), $dom#5616 <= $Fluid#5618
do
  for $c#5621 : int3d(Fluid_columns, $dom#5616) in $dom#5616 do
    $Fluid#5618[$c#5621].rho_t = double(double(0))
    $Fluid#5618[$c#5621].rhoVelocity_t = double[3](array(double(0), double(0), double(0)))
    $Fluid#5618[$c#5621].rhoEnergy_t = double(double(0))
    $Fluid#5618[$c#5621].rhoEnthalpy = ($Fluid#5618[$c#5621].rhoEnergy+$Fluid#5618[$c#5621].pressure)
  end
end
task Particles_InitializeTimeDerivatives($dom#5651 : region#1200(ispace#1200(int1d), particles_columns), $particles#5653 : region#1201(ispace#1201(int1d), particles_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($particles#5653.position_t), writes($particles#5653.position_t), reads($particles#5653.temperature_t), writes($particles#5653.temperature_t), reads($particles#5653.velocity_t), writes($particles#5653.velocity_t), reads($particles#5653.__valid), $dom#5651 <= $particles#5653
do
  for $p#5656 : int1d(particles_columns, $dom#5651) in $dom#5651 do
    if $particles#5653[$p#5656].__valid then
      $particles#5653[$p#5656].position_t = double[3](array(double(0), double(0), double(0)))
      $particles#5653[$p#5656].velocity_t = double[3](array(double(0), double(0), double(0)))
      $particles#5653[$p#5656].temperature_t = double(int32(0))
    else
    end
  end
end
task Flow_UpdateGhostVelocityGradientStep1($dom#5692 : region#1211(ispace#1211(int3d), Fluid_columns), $Fluid#5694 : region#1212(ispace#1212(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#5694.velocityGradientX), reads($Fluid#5694.velocityGradientXBoundary), writes($Fluid#5694.velocityGradientXBoundary), reads($Fluid#5694.velocityGradientY), reads($Fluid#5694.velocityGradientYBoundary), writes($Fluid#5694.velocityGradientYBoundary), reads($Fluid#5694.velocityGradientZ), reads($Fluid#5694.velocityGradientZBoundary), writes($Fluid#5694.velocityGradientZBoundary), $dom#5692 <= $Fluid#5694
do
  for $c#5823 : int3d(Fluid_columns, $dom#5692) in $dom#5692 do
    if (max(int32((uint64(int32(0))-int3d($c#5823).x)), int32(0))>int32(0)) then
      var $c_bnd#5824 : int3d(Fluid_columns, $dom#5692) = $c#5823
      var $c_int#5825 : int3d = (($c#5823+{1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#5826 : double[3] = array(double(1), double(1), double(1))
      $Fluid#5694[$c_bnd#5824].velocityGradientXBoundary = vv_mul_double_3($sign#5826, $Fluid#5694[$c_int#5825].velocityGradientX)
      $Fluid#5694[$c_bnd#5824].velocityGradientYBoundary = vv_mul_double_3($sign#5826, $Fluid#5694[$c_int#5825].velocityGradientY)
      $Fluid#5694[$c_bnd#5824].velocityGradientZBoundary = vv_mul_double_3($sign#5826, $Fluid#5694[$c_int#5825].velocityGradientZ)
    else
    end
    if (max(int32((int3d($c#5823).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0)) then
      var $c_bnd#5827 : int3d(Fluid_columns, $dom#5692) = $c#5823
      var $c_int#5828 : int3d = (($c#5823+{-1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#5829 : double[3] = array(double(1), double(1), double(1))
      $Fluid#5694[$c_bnd#5827].velocityGradientXBoundary = vv_mul_double_3($sign#5829, $Fluid#5694[$c_int#5828].velocityGradientX)
      $Fluid#5694[$c_bnd#5827].velocityGradientYBoundary = vv_mul_double_3($sign#5829, $Fluid#5694[$c_int#5828].velocityGradientY)
      $Fluid#5694[$c_bnd#5827].velocityGradientZBoundary = vv_mul_double_3($sign#5829, $Fluid#5694[$c_int#5828].velocityGradientZ)
    else
    end
    if (max(int32((uint64(int32(1))-int3d($c#5823).y)), int32(0))>int32(0)) then
      var $c_bnd#5830 : int3d(Fluid_columns, $dom#5692) = $c#5823
      var $c_int#5831 : int3d = (($c#5823+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#5832 : double[3] = array(double(-1), double(-1), double(-1))
      $Fluid#5694[$c_bnd#5830].velocityGradientXBoundary = vv_mul_double_3($sign#5832, $Fluid#5694[$c_int#5831].velocityGradientX)
      $Fluid#5694[$c_bnd#5830].velocityGradientYBoundary = vv_mul_double_3($sign#5832, $Fluid#5694[$c_int#5831].velocityGradientY)
      $Fluid#5694[$c_bnd#5830].velocityGradientZBoundary = vv_mul_double_3($sign#5832, $Fluid#5694[$c_int#5831].velocityGradientZ)
    else
    end
    if (max(int32((int3d($c#5823).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0)) then
      var $c_bnd#5833 : int3d(Fluid_columns, $dom#5692) = $c#5823
      var $c_int#5834 : int3d = (($c#5823+{0, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#5835 : double[3] = array(double(-1), double(-1), double(-1))
      $Fluid#5694[$c_bnd#5833].velocityGradientXBoundary = vv_mul_double_3($sign#5835, $Fluid#5694[$c_int#5834].velocityGradientX)
      $Fluid#5694[$c_bnd#5833].velocityGradientYBoundary = vv_mul_double_3($sign#5835, $Fluid#5694[$c_int#5834].velocityGradientY)
      $Fluid#5694[$c_bnd#5833].velocityGradientZBoundary = vv_mul_double_3($sign#5835, $Fluid#5694[$c_int#5834].velocityGradientZ)
    else
    end
    if (max(int32((uint64(int32(1))-int3d($c#5823).z)), int32(0))>int32(0)) then
      var $c_bnd#5836 : int3d(Fluid_columns, $dom#5692) = $c#5823
      var $c_int#5837 : int3d = (($c#5823+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#5838 : double[3] = array(double(1), double(1), double(-1))
      $Fluid#5694[$c_bnd#5836].velocityGradientXBoundary = vv_mul_double_3($sign#5838, $Fluid#5694[$c_int#5837].velocityGradientX)
      $Fluid#5694[$c_bnd#5836].velocityGradientYBoundary = vv_mul_double_3($sign#5838, $Fluid#5694[$c_int#5837].velocityGradientY)
      $Fluid#5694[$c_bnd#5836].velocityGradientZBoundary = vv_mul_double_3($sign#5838, $Fluid#5694[$c_int#5837].velocityGradientZ)
    else
    end
    if (max(int32((int3d($c#5823).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0)) then
      var $c_bnd#5839 : int3d(Fluid_columns, $dom#5692) = $c#5823
      var $c_int#5840 : int3d = (($c#5823+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
      var $sign#5841 : double[3] = array(double(1), double(1), double(-1))
      $Fluid#5694[$c_bnd#5839].velocityGradientXBoundary = vv_mul_double_3($sign#5841, $Fluid#5694[$c_int#5840].velocityGradientX)
      $Fluid#5694[$c_bnd#5839].velocityGradientYBoundary = vv_mul_double_3($sign#5841, $Fluid#5694[$c_int#5840].velocityGradientY)
      $Fluid#5694[$c_bnd#5839].velocityGradientZBoundary = vv_mul_double_3($sign#5841, $Fluid#5694[$c_int#5840].velocityGradientZ)
    else
    end
  end
end
task Flow_UpdateGhostVelocityGradientStep2($dom#6111 : region#1306(ispace#1306(int3d), Fluid_columns), $Fluid#6113 : region#1307(ispace#1307(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#6113.velocityGradientX), writes($Fluid#6113.velocityGradientX), reads($Fluid#6113.velocityGradientXBoundary), reads($Fluid#6113.velocityGradientY), writes($Fluid#6113.velocityGradientY), reads($Fluid#6113.velocityGradientYBoundary), reads($Fluid#6113.velocityGradientZ), writes($Fluid#6113.velocityGradientZ), reads($Fluid#6113.velocityGradientZBoundary), $dom#6111 <= $Fluid#6113
do
  for $c#6116 : int3d(Fluid_columns, $dom#6111) in $dom#6111 do
    if ((((((max(int32((uint64(int32(0))-int3d($c#6116).x)), int32(0))>int32(0)) or (max(int32((int3d($c#6116).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#6116).y)), int32(0))>int32(0))) or (max(int32((int3d($c#6116).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#6116).z)), int32(0))>int32(0))) or (max(int32((int3d($c#6116).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0))) then
      $Fluid#6113[$c#6116].velocityGradientX = $Fluid#6113[$c#6116].velocityGradientXBoundary
      $Fluid#6113[$c#6116].velocityGradientY = $Fluid#6113[$c#6116].velocityGradientYBoundary
      $Fluid#6113[$c#6116].velocityGradientZ = $Fluid#6113[$c#6116].velocityGradientZBoundary
    else
    end
  end
end
task CenteredInviscidFlux($c_l#6145 : int3d, $c_r#6146 : int3d, $Fluid#6148 : region#1319(ispace#1319(int3d), Fluid_columns)) : double[5]
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#6148.pressure), reads($Fluid#6148.rho), reads($Fluid#6148.rhoEnthalpy), reads($Fluid#6148.rhoVelocity), reads($Fluid#6148.velocity)
do
  var $rhoFactorDiagonal#6201 : double = double(double(0))
  var $rhoVelocityFactorDiagonal#6202 : double[3] = double[3](array(double(0), double(0), double(0)))
  var $rhoEnergyFactorDiagonal#6203 : double = double(double(0))
  var $fpdiag#6204 : double = double(double(0))
  $rhoFactorDiagonal#6201 = (double(0.5)*(($Fluid#6148[$c_l#6145].rho*$Fluid#6148[$c_l#6145].velocity[int32(0)])+($Fluid#6148[$c_r#6146].rho*$Fluid#6148[$c_r#6146].velocity[int32(0)])))
  $rhoVelocityFactorDiagonal#6202 = vs_mul_double_3(vv_add_double_3(vs_mul_double_3($Fluid#6148[$c_l#6145].rhoVelocity, $Fluid#6148[$c_l#6145].velocity[int32(0)]), vs_mul_double_3($Fluid#6148[$c_r#6146].rhoVelocity, $Fluid#6148[$c_r#6146].velocity[int32(0)])), double(0.5))
  $rhoEnergyFactorDiagonal#6203 = (double(0.5)*(($Fluid#6148[$c_l#6145].rhoEnthalpy*$Fluid#6148[$c_l#6145].velocity[int32(0)])+($Fluid#6148[$c_r#6146].rhoEnthalpy*$Fluid#6148[$c_r#6146].velocity[int32(0)])))
  $fpdiag#6204 += (double(0.5)*($Fluid#6148[$c_l#6145].pressure+$Fluid#6148[$c_r#6146].pressure))
  var $rhoFactorSkew#6205 : double = double(double(0))
  var $rhoVelocityFactorSkew#6206 : double[3] = double[3](array(double(0), double(0), double(0)))
  var $rhoEnergyFactorSkew#6207 : double = double(double(0))
  var $tmp#6208 : double = double(double(0))
  $tmp#6208 = (double(0.5)*$Fluid#6148[$c_r#6146].velocity[int32(0)])
  $rhoFactorSkew#6205 += ($Fluid#6148[$c_l#6145].rho*$tmp#6208)
  var $tmp#6209 : double[3] = vs_mul_double_3($Fluid#6148[$c_l#6145].rhoVelocity, $tmp#6208)
  var $v#6210 : double[3] = $rhoVelocityFactorSkew#6206
  $v#6210[0] += $tmp#6209[0]
  $v#6210[1] += $tmp#6209[1]
  $v#6210[2] += $tmp#6209[2]
  $rhoVelocityFactorSkew#6206 = $v#6210
  $rhoEnergyFactorSkew#6207 += ($Fluid#6148[$c_l#6145].rhoEnthalpy*$tmp#6208)
  $tmp#6208 = (double(0.5)*$Fluid#6148[$c_l#6145].velocity[int32(0)])
  $rhoFactorSkew#6205 += ($Fluid#6148[$c_r#6146].rho*$tmp#6208)
  var $tmp#6211 : double[3] = vs_mul_double_3($Fluid#6148[$c_r#6146].rhoVelocity, $tmp#6208)
  var $v#6212 : double[3] = $rhoVelocityFactorSkew#6206
  $v#6212[0] += $tmp#6211[0]
  $v#6212[1] += $tmp#6211[1]
  $v#6212[2] += $tmp#6211[2]
  $rhoVelocityFactorSkew#6206 = $v#6212
  $rhoEnergyFactorSkew#6207 += ($Fluid#6148[$c_r#6146].rhoEnthalpy*$tmp#6208)
  var $s#6213 : double = double(0.5)
  var $rhoFlux_temp#6214 : double = (($s#6213*$rhoFactorDiagonal#6201)+((double(int32(1))-$s#6213)*$rhoFactorSkew#6205))
  var $rhoVelocityFlux_temp#6215 : double[3] = vv_add_double_3(vs_mul_double_3($rhoVelocityFactorDiagonal#6202, $s#6213), vs_mul_double_3($rhoVelocityFactorSkew#6206, (double(int32(1))-$s#6213)))
  var $rhoEnergyFlux_temp#6216 : double = (($s#6213*$rhoEnergyFactorDiagonal#6203)+((double(int32(1))-$s#6213)*$rhoEnergyFactorSkew#6207))
  $rhoVelocityFlux_temp#6215[int32(0)] += $fpdiag#6204
  return array($rhoFlux_temp#6214, $rhoVelocityFlux_temp#6215[int32(0)], $rhoVelocityFlux_temp#6215[int32(1)], $rhoVelocityFlux_temp#6215[int32(2)], $rhoEnergyFlux_temp#6216)
end
task CenteredInviscidFlux_($c_l#6565 : int3d, $c_r#6566 : int3d, $Fluid#6568 : region#1435(ispace#1435(int3d), Fluid_columns)) : double[5]
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#6568.pressure), reads($Fluid#6568.rho), reads($Fluid#6568.rhoEnthalpy), reads($Fluid#6568.rhoVelocity), reads($Fluid#6568.velocity)
do
  var $rhoFactorDiagonal#6621 : double = double(double(0))
  var $rhoVelocityFactorDiagonal#6622 : double[3] = double[3](array(double(0), double(0), double(0)))
  var $rhoEnergyFactorDiagonal#6623 : double = double(double(0))
  var $fpdiag#6624 : double = double(double(0))
  $rhoFactorDiagonal#6621 = (double(0.5)*(($Fluid#6568[$c_l#6565].rho*$Fluid#6568[$c_l#6565].velocity[int32(1)])+($Fluid#6568[$c_r#6566].rho*$Fluid#6568[$c_r#6566].velocity[int32(1)])))
  $rhoVelocityFactorDiagonal#6622 = vs_mul_double_3(vv_add_double_3(vs_mul_double_3($Fluid#6568[$c_l#6565].rhoVelocity, $Fluid#6568[$c_l#6565].velocity[int32(1)]), vs_mul_double_3($Fluid#6568[$c_r#6566].rhoVelocity, $Fluid#6568[$c_r#6566].velocity[int32(1)])), double(0.5))
  $rhoEnergyFactorDiagonal#6623 = (double(0.5)*(($Fluid#6568[$c_l#6565].rhoEnthalpy*$Fluid#6568[$c_l#6565].velocity[int32(1)])+($Fluid#6568[$c_r#6566].rhoEnthalpy*$Fluid#6568[$c_r#6566].velocity[int32(1)])))
  $fpdiag#6624 += (double(0.5)*($Fluid#6568[$c_l#6565].pressure+$Fluid#6568[$c_r#6566].pressure))
  var $rhoFactorSkew#6625 : double = double(double(0))
  var $rhoVelocityFactorSkew#6626 : double[3] = double[3](array(double(0), double(0), double(0)))
  var $rhoEnergyFactorSkew#6627 : double = double(double(0))
  var $tmp#6628 : double = double(double(0))
  $tmp#6628 = (double(0.5)*$Fluid#6568[$c_r#6566].velocity[int32(1)])
  $rhoFactorSkew#6625 += ($Fluid#6568[$c_l#6565].rho*$tmp#6628)
  var $tmp#6629 : double[3] = vs_mul_double_3($Fluid#6568[$c_l#6565].rhoVelocity, $tmp#6628)
  var $v#6630 : double[3] = $rhoVelocityFactorSkew#6626
  $v#6630[0] += $tmp#6629[0]
  $v#6630[1] += $tmp#6629[1]
  $v#6630[2] += $tmp#6629[2]
  $rhoVelocityFactorSkew#6626 = $v#6630
  $rhoEnergyFactorSkew#6627 += ($Fluid#6568[$c_l#6565].rhoEnthalpy*$tmp#6628)
  $tmp#6628 = (double(0.5)*$Fluid#6568[$c_l#6565].velocity[int32(1)])
  $rhoFactorSkew#6625 += ($Fluid#6568[$c_r#6566].rho*$tmp#6628)
  var $tmp#6631 : double[3] = vs_mul_double_3($Fluid#6568[$c_r#6566].rhoVelocity, $tmp#6628)
  var $v#6632 : double[3] = $rhoVelocityFactorSkew#6626
  $v#6632[0] += $tmp#6631[0]
  $v#6632[1] += $tmp#6631[1]
  $v#6632[2] += $tmp#6631[2]
  $rhoVelocityFactorSkew#6626 = $v#6632
  $rhoEnergyFactorSkew#6627 += ($Fluid#6568[$c_r#6566].rhoEnthalpy*$tmp#6628)
  var $s#6633 : double = double(0.5)
  var $rhoFlux_temp#6634 : double = (($s#6633*$rhoFactorDiagonal#6621)+((double(int32(1))-$s#6633)*$rhoFactorSkew#6625))
  var $rhoVelocityFlux_temp#6635 : double[3] = vv_add_double_3(vs_mul_double_3($rhoVelocityFactorDiagonal#6622, $s#6633), vs_mul_double_3($rhoVelocityFactorSkew#6626, (double(int32(1))-$s#6633)))
  var $rhoEnergyFlux_temp#6636 : double = (($s#6633*$rhoEnergyFactorDiagonal#6623)+((double(int32(1))-$s#6633)*$rhoEnergyFactorSkew#6627))
  $rhoVelocityFlux_temp#6635[int32(1)] += $fpdiag#6624
  return array($rhoFlux_temp#6634, $rhoVelocityFlux_temp#6635[int32(0)], $rhoVelocityFlux_temp#6635[int32(1)], $rhoVelocityFlux_temp#6635[int32(2)], $rhoEnergyFlux_temp#6636)
end
task CenteredInviscidFlux__($c_l#6985 : int3d, $c_r#6986 : int3d, $Fluid#6988 : region#1551(ispace#1551(int3d), Fluid_columns)) : double[5]
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#6988.pressure), reads($Fluid#6988.rho), reads($Fluid#6988.rhoEnthalpy), reads($Fluid#6988.rhoVelocity), reads($Fluid#6988.velocity)
do
  var $rhoFactorDiagonal#7041 : double = double(double(0))
  var $rhoVelocityFactorDiagonal#7042 : double[3] = double[3](array(double(0), double(0), double(0)))
  var $rhoEnergyFactorDiagonal#7043 : double = double(double(0))
  var $fpdiag#7044 : double = double(double(0))
  $rhoFactorDiagonal#7041 = (double(0.5)*(($Fluid#6988[$c_l#6985].rho*$Fluid#6988[$c_l#6985].velocity[int32(2)])+($Fluid#6988[$c_r#6986].rho*$Fluid#6988[$c_r#6986].velocity[int32(2)])))
  $rhoVelocityFactorDiagonal#7042 = vs_mul_double_3(vv_add_double_3(vs_mul_double_3($Fluid#6988[$c_l#6985].rhoVelocity, $Fluid#6988[$c_l#6985].velocity[int32(2)]), vs_mul_double_3($Fluid#6988[$c_r#6986].rhoVelocity, $Fluid#6988[$c_r#6986].velocity[int32(2)])), double(0.5))
  $rhoEnergyFactorDiagonal#7043 = (double(0.5)*(($Fluid#6988[$c_l#6985].rhoEnthalpy*$Fluid#6988[$c_l#6985].velocity[int32(2)])+($Fluid#6988[$c_r#6986].rhoEnthalpy*$Fluid#6988[$c_r#6986].velocity[int32(2)])))
  $fpdiag#7044 += (double(0.5)*($Fluid#6988[$c_l#6985].pressure+$Fluid#6988[$c_r#6986].pressure))
  var $rhoFactorSkew#7045 : double = double(double(0))
  var $rhoVelocityFactorSkew#7046 : double[3] = double[3](array(double(0), double(0), double(0)))
  var $rhoEnergyFactorSkew#7047 : double = double(double(0))
  var $tmp#7048 : double = double(double(0))
  $tmp#7048 = (double(0.5)*$Fluid#6988[$c_r#6986].velocity[int32(2)])
  $rhoFactorSkew#7045 += ($Fluid#6988[$c_l#6985].rho*$tmp#7048)
  var $tmp#7049 : double[3] = vs_mul_double_3($Fluid#6988[$c_l#6985].rhoVelocity, $tmp#7048)
  var $v#7050 : double[3] = $rhoVelocityFactorSkew#7046
  $v#7050[0] += $tmp#7049[0]
  $v#7050[1] += $tmp#7049[1]
  $v#7050[2] += $tmp#7049[2]
  $rhoVelocityFactorSkew#7046 = $v#7050
  $rhoEnergyFactorSkew#7047 += ($Fluid#6988[$c_l#6985].rhoEnthalpy*$tmp#7048)
  $tmp#7048 = (double(0.5)*$Fluid#6988[$c_l#6985].velocity[int32(2)])
  $rhoFactorSkew#7045 += ($Fluid#6988[$c_r#6986].rho*$tmp#7048)
  var $tmp#7051 : double[3] = vs_mul_double_3($Fluid#6988[$c_r#6986].rhoVelocity, $tmp#7048)
  var $v#7052 : double[3] = $rhoVelocityFactorSkew#7046
  $v#7052[0] += $tmp#7051[0]
  $v#7052[1] += $tmp#7051[1]
  $v#7052[2] += $tmp#7051[2]
  $rhoVelocityFactorSkew#7046 = $v#7052
  $rhoEnergyFactorSkew#7047 += ($Fluid#6988[$c_r#6986].rhoEnthalpy*$tmp#7048)
  var $s#7053 : double = double(0.5)
  var $rhoFlux_temp#7054 : double = (($s#7053*$rhoFactorDiagonal#7041)+((double(int32(1))-$s#7053)*$rhoFactorSkew#7045))
  var $rhoVelocityFlux_temp#7055 : double[3] = vv_add_double_3(vs_mul_double_3($rhoVelocityFactorDiagonal#7042, $s#7053), vs_mul_double_3($rhoVelocityFactorSkew#7046, (double(int32(1))-$s#7053)))
  var $rhoEnergyFlux_temp#7056 : double = (($s#7053*$rhoEnergyFactorDiagonal#7043)+((double(int32(1))-$s#7053)*$rhoEnergyFactorSkew#7047))
  $rhoVelocityFlux_temp#7055[int32(2)] += $fpdiag#7044
  return array($rhoFlux_temp#7054, $rhoVelocityFlux_temp#7055[int32(0)], $rhoVelocityFlux_temp#7055[int32(1)], $rhoVelocityFlux_temp#7055[int32(2)], $rhoEnergyFlux_temp#7056)
end
task Flow_AddGetFlux($dom#6140 : region#1317(ispace#1317(int3d), Fluid_columns), $Fluid#6142 : region#1318(ispace#1318(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#6142.pressure), reads($Fluid#6142.rho), reads($Fluid#6142.rhoEnergyFluxX), writes($Fluid#6142.rhoEnergyFluxX), reads($Fluid#6142.rhoEnergyFluxX), writes($Fluid#6142.rhoEnergyFluxX), reads($Fluid#6142.rhoEnergyFluxY), writes($Fluid#6142.rhoEnergyFluxY), reads($Fluid#6142.rhoEnergyFluxY), writes($Fluid#6142.rhoEnergyFluxY), reads($Fluid#6142.rhoEnergyFluxZ), writes($Fluid#6142.rhoEnergyFluxZ), reads($Fluid#6142.rhoEnergyFluxZ), writes($Fluid#6142.rhoEnergyFluxZ), reads($Fluid#6142.rhoEnthalpy), reads($Fluid#6142.rhoFluxX), writes($Fluid#6142.rhoFluxX), reads($Fluid#6142.rhoFluxY), writes($Fluid#6142.rhoFluxY), reads($Fluid#6142.rhoFluxZ), writes($Fluid#6142.rhoFluxZ), reads($Fluid#6142.rhoVelocity), reads($Fluid#6142.rhoVelocityFluxX), writes($Fluid#6142.rhoVelocityFluxX), reads($Fluid#6142.rhoVelocityFluxX), writes($Fluid#6142.rhoVelocityFluxX), reads($Fluid#6142.rhoVelocityFluxY), writes($Fluid#6142.rhoVelocityFluxY), reads($Fluid#6142.rhoVelocityFluxY), writes($Fluid#6142.rhoVelocityFluxY), reads($Fluid#6142.rhoVelocityFluxZ), writes($Fluid#6142.rhoVelocityFluxZ), reads($Fluid#6142.rhoVelocityFluxZ), writes($Fluid#6142.rhoVelocityFluxZ), reads($Fluid#6142.temperature), reads($Fluid#6142.velocity), reads($Fluid#6142.velocityGradientX), reads($Fluid#6142.velocityGradientY), reads($Fluid#6142.velocityGradientZ), $dom#6140 <= $Fluid#6142
do
  for $c#7529 : int3d(Fluid_columns, $dom#6140) in $dom#6140 do
    if ((not ((((((max(int32((uint64(int32(0))-int3d($c#7529).x)), int32(0))>int32(0)) or (max(int32((int3d($c#7529).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#7529).y)), int32(0))>int32(0))) or (max(int32((int3d($c#7529).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#7529).z)), int32(0))>int32(0))) or (max(int32((int3d($c#7529).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0)))) or (max(int32((uint64(int32(0))-int3d($c#7529).x)), int32(0))==int32(1))) then
      var $7712 : double[5]
      do
        var $7710 : int3d = int3d($c#7529)
        var $7711 : int3d = (($c#7529+{1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
        var $rhoFactorDiagonal#6201 : double = double(double(0))
        var $rhoVelocityFactorDiagonal#6202 : double[3] = double[3](array(double(0), double(0), double(0)))
        var $rhoEnergyFactorDiagonal#6203 : double = double(double(0))
        var $fpdiag#6204 : double = double(double(0))
        $rhoFactorDiagonal#6201 = (double(0.5)*(($Fluid#6142[$7710].rho*$Fluid#6142[$7710].velocity[int32(0)])+($Fluid#6142[$7711].rho*$Fluid#6142[$7711].velocity[int32(0)])))
        $rhoVelocityFactorDiagonal#6202 = vs_mul_double_3(vv_add_double_3(vs_mul_double_3($Fluid#6142[$7710].rhoVelocity, $Fluid#6142[$7710].velocity[int32(0)]), vs_mul_double_3($Fluid#6142[$7711].rhoVelocity, $Fluid#6142[$7711].velocity[int32(0)])), double(0.5))
        $rhoEnergyFactorDiagonal#6203 = (double(0.5)*(($Fluid#6142[$7710].rhoEnthalpy*$Fluid#6142[$7710].velocity[int32(0)])+($Fluid#6142[$7711].rhoEnthalpy*$Fluid#6142[$7711].velocity[int32(0)])))
        $fpdiag#6204 += (double(0.5)*($Fluid#6142[$7710].pressure+$Fluid#6142[$7711].pressure))
        var $rhoFactorSkew#6205 : double = double(double(0))
        var $rhoVelocityFactorSkew#6206 : double[3] = double[3](array(double(0), double(0), double(0)))
        var $rhoEnergyFactorSkew#6207 : double = double(double(0))
        var $tmp#6208 : double = double(double(0))
        $tmp#6208 = (double(0.5)*$Fluid#6142[$7711].velocity[int32(0)])
        $rhoFactorSkew#6205 += ($Fluid#6142[$7710].rho*$tmp#6208)
        var $tmp#6209 : double[3] = vs_mul_double_3($Fluid#6142[$7710].rhoVelocity, $tmp#6208)
        var $v#6210 : double[3] = $rhoVelocityFactorSkew#6206
        $v#6210[0] += $tmp#6209[0]
        $v#6210[1] += $tmp#6209[1]
        $v#6210[2] += $tmp#6209[2]
        $rhoVelocityFactorSkew#6206 = $v#6210
        $rhoEnergyFactorSkew#6207 += ($Fluid#6142[$7710].rhoEnthalpy*$tmp#6208)
        $tmp#6208 = (double(0.5)*$Fluid#6142[$7710].velocity[int32(0)])
        $rhoFactorSkew#6205 += ($Fluid#6142[$7711].rho*$tmp#6208)
        var $tmp#6211 : double[3] = vs_mul_double_3($Fluid#6142[$7711].rhoVelocity, $tmp#6208)
        var $v#6212 : double[3] = $rhoVelocityFactorSkew#6206
        $v#6212[0] += $tmp#6211[0]
        $v#6212[1] += $tmp#6211[1]
        $v#6212[2] += $tmp#6211[2]
        $rhoVelocityFactorSkew#6206 = $v#6212
        $rhoEnergyFactorSkew#6207 += ($Fluid#6142[$7711].rhoEnthalpy*$tmp#6208)
        var $s#6213 : double = double(0.5)
        var $rhoFlux_temp#6214 : double = (($s#6213*$rhoFactorDiagonal#6201)+((double(int32(1))-$s#6213)*$rhoFactorSkew#6205))
        var $rhoVelocityFlux_temp#6215 : double[3] = vv_add_double_3(vs_mul_double_3($rhoVelocityFactorDiagonal#6202, $s#6213), vs_mul_double_3($rhoVelocityFactorSkew#6206, (double(int32(1))-$s#6213)))
        var $rhoEnergyFlux_temp#6216 : double = (($s#6213*$rhoEnergyFactorDiagonal#6203)+((double(int32(1))-$s#6213)*$rhoEnergyFactorSkew#6207))
        $rhoVelocityFlux_temp#6215[int32(0)] += $fpdiag#6204
        $7712 = array($rhoFlux_temp#6214, $rhoVelocityFlux_temp#6215[int32(0)], $rhoVelocityFlux_temp#6215[int32(1)], $rhoVelocityFlux_temp#6215[int32(2)], $rhoEnergyFlux_temp#6216)
      end
      var $dummy#7713 : int32 = 0
      var $flux#7530 : double[5] = $7712
      $Fluid#6142[$c#7529].rhoFluxX = $flux#7530[int32(0)]
      $Fluid#6142[$c#7529].rhoVelocityFluxX = array($flux#7530[int32(1)], $flux#7530[int32(2)], $flux#7530[int32(3)])
      $Fluid#6142[$c#7529].rhoEnergyFluxX = $flux#7530[int32(4)]
      var $7715 : double
      do
        var $7714 : double = $Fluid#6142[(($c#7529+{1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        var $viscosity#5391 : double = double(double(0))
        if (int32(3003)==int32(3001)) then
          $viscosity#5391 = double(0.2454194)
        else
          if (int32(3003)==int32(3002)) then
            $viscosity#5391 = (double(0.9683056)*pow(($7714/double(0.9913615)), double(0.75)))
          else
            if (int32(3003)==int32(3003)) then
              $viscosity#5391 = ((double(0.5815515)*pow(($7714/double(0.6935851)), (double(3)/double(2))))*((double(0.6935851)+double(0.2772065))/($7714+double(0.2772065))))
            else
              std.assert(false, "(Liszt assertion)")
            end
          end
        end
        $7715 = $viscosity#5391
      end
      var $dummy#7716 : int32 = 0
      var $7718 : double
      do
        var $7717 : double = $Fluid#6142[$c#7529].temperature
        var $viscosity#5391 : double = double(double(0))
        if (int32(3003)==int32(3001)) then
          $viscosity#5391 = double(0.2454194)
        else
          if (int32(3003)==int32(3002)) then
            $viscosity#5391 = (double(0.9683056)*pow(($7717/double(0.9913615)), double(0.75)))
          else
            if (int32(3003)==int32(3003)) then
              $viscosity#5391 = ((double(0.5815515)*pow(($7717/double(0.6935851)), (double(3)/double(2))))*((double(0.6935851)+double(0.2772065))/($7717+double(0.2772065))))
            else
              std.assert(false, "(Liszt assertion)")
            end
          end
        end
        $7718 = $viscosity#5391
      end
      var $dummy#7719 : int32 = 0
      var $muFace#7531 : double = (double(0.5)*($7718+$7715))
      var $velocityFace#7532 : double[3] = double[3](array(double(0), double(0), double(0)))
      var $velocityX_YFace#7533 : double = double(double(0))
      var $velocityX_ZFace#7534 : double = double(double(0))
      var $velocityY_YFace#7535 : double = double(double(0))
      var $velocityZ_ZFace#7536 : double = double(double(0))
      $velocityFace#7532 = vs_mul_double_3(vv_add_double_3($Fluid#6142[$c#7529].velocity, $Fluid#6142[(($c#7529+{1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity), double(0.5))
      $velocityX_YFace#7533 = (double(0.5)*($Fluid#6142[$c#7529].velocityGradientY[int32(0)]+$Fluid#6142[(($c#7529+{1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocityGradientY[int32(0)]))
      $velocityX_ZFace#7534 = (double(0.5)*($Fluid#6142[$c#7529].velocityGradientZ[int32(0)]+$Fluid#6142[(($c#7529+{1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocityGradientZ[int32(0)]))
      $velocityY_YFace#7535 = (double(0.5)*($Fluid#6142[$c#7529].velocityGradientY[int32(1)]+$Fluid#6142[(($c#7529+{1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocityGradientY[int32(1)]))
      $velocityZ_ZFace#7536 = (double(0.5)*($Fluid#6142[$c#7529].velocityGradientZ[int32(2)]+$Fluid#6142[(($c#7529+{1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocityGradientZ[int32(2)]))
      var $velocityX_XFace#7537 : double = double(double(0))
      var $velocityY_XFace#7538 : double = double(double(0))
      var $velocityZ_XFace#7539 : double = double(double(0))
      var $temperature_XFace#7540 : double = double(double(0))
      $velocityX_XFace#7537 = (double(0.5)*($Fluid#6142[(($c#7529+{1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity[int32(0)]-$Fluid#6142[$c#7529].velocity[int32(0)]))
      $velocityY_XFace#7538 = (double(0.5)*($Fluid#6142[(($c#7529+{1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity[int32(1)]-$Fluid#6142[$c#7529].velocity[int32(1)]))
      $velocityZ_XFace#7539 = (double(0.5)*($Fluid#6142[(($c#7529+{1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity[int32(2)]-$Fluid#6142[$c#7529].velocity[int32(2)]))
      $temperature_XFace#7540 = (double(0.5)*($Fluid#6142[(($c#7529+{1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature-$Fluid#6142[$c#7529].temperature))
      $velocityX_XFace#7537 *= (1/(double(0.0015841600529101)*double(0.5)))
      $velocityY_XFace#7538 *= (1/(double(0.0015841600529101)*double(0.5)))
      $velocityZ_XFace#7539 *= (1/(double(0.0015841600529101)*double(0.5)))
      $temperature_XFace#7540 *= (1/(double(0.0015841600529101)*double(0.5)))
      var $sigmaXX#7541 : double = (($muFace#7531*(((double(4)*$velocityX_XFace#7537)-(double(2)*$velocityY_YFace#7535))-(double(2)*$velocityZ_ZFace#7536)))/double(3))
      var $sigmaYX#7542 : double = ($muFace#7531*($velocityY_XFace#7538+$velocityX_YFace#7533))
      var $sigmaZX#7543 : double = ($muFace#7531*($velocityZ_XFace#7539+$velocityX_ZFace#7534))
      var $usigma#7544 : double = ((($velocityFace#7532[int32(0)]*$sigmaXX#7541)+($velocityFace#7532[int32(1)]*$sigmaYX#7542))+($velocityFace#7532[int32(2)]*$sigmaZX#7543))
      var $cp#7545 : double = ((double(0.8198141)*double(0.7272462))/(double(0.8198141)-double(1)))
      var $heatFlux#7546 : double = ((-(($cp#7545*$muFace#7531)/double(0.1565651)))*$temperature_XFace#7540)
      $Fluid#6142[$c#7529].rhoVelocityFluxX[int32(0)] += (-$sigmaXX#7541)
      $Fluid#6142[$c#7529].rhoVelocityFluxX[int32(1)] += (-$sigmaYX#7542)
      $Fluid#6142[$c#7529].rhoVelocityFluxX[int32(2)] += (-$sigmaZX#7543)
      $Fluid#6142[$c#7529].rhoEnergyFluxX += (-($usigma#7544-$heatFlux#7546))
    else
    end
    if ((not ((((((max(int32((uint64(int32(0))-int3d($c#7529).x)), int32(0))>int32(0)) or (max(int32((int3d($c#7529).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#7529).y)), int32(0))>int32(0))) or (max(int32((int3d($c#7529).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#7529).z)), int32(0))>int32(0))) or (max(int32((int3d($c#7529).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0)))) or (max(int32((uint64(int32(1))-int3d($c#7529).y)), int32(0))==int32(1))) then
      var $7722 : double[5]
      do
        var $7720 : int3d = int3d($c#7529)
        var $7721 : int3d = (($c#7529+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
        var $rhoFactorDiagonal#6621 : double = double(double(0))
        var $rhoVelocityFactorDiagonal#6622 : double[3] = double[3](array(double(0), double(0), double(0)))
        var $rhoEnergyFactorDiagonal#6623 : double = double(double(0))
        var $fpdiag#6624 : double = double(double(0))
        $rhoFactorDiagonal#6621 = (double(0.5)*(($Fluid#6142[$7720].rho*$Fluid#6142[$7720].velocity[int32(1)])+($Fluid#6142[$7721].rho*$Fluid#6142[$7721].velocity[int32(1)])))
        $rhoVelocityFactorDiagonal#6622 = vs_mul_double_3(vv_add_double_3(vs_mul_double_3($Fluid#6142[$7720].rhoVelocity, $Fluid#6142[$7720].velocity[int32(1)]), vs_mul_double_3($Fluid#6142[$7721].rhoVelocity, $Fluid#6142[$7721].velocity[int32(1)])), double(0.5))
        $rhoEnergyFactorDiagonal#6623 = (double(0.5)*(($Fluid#6142[$7720].rhoEnthalpy*$Fluid#6142[$7720].velocity[int32(1)])+($Fluid#6142[$7721].rhoEnthalpy*$Fluid#6142[$7721].velocity[int32(1)])))
        $fpdiag#6624 += (double(0.5)*($Fluid#6142[$7720].pressure+$Fluid#6142[$7721].pressure))
        var $rhoFactorSkew#6625 : double = double(double(0))
        var $rhoVelocityFactorSkew#6626 : double[3] = double[3](array(double(0), double(0), double(0)))
        var $rhoEnergyFactorSkew#6627 : double = double(double(0))
        var $tmp#6628 : double = double(double(0))
        $tmp#6628 = (double(0.5)*$Fluid#6142[$7721].velocity[int32(1)])
        $rhoFactorSkew#6625 += ($Fluid#6142[$7720].rho*$tmp#6628)
        var $tmp#6629 : double[3] = vs_mul_double_3($Fluid#6142[$7720].rhoVelocity, $tmp#6628)
        var $v#6630 : double[3] = $rhoVelocityFactorSkew#6626
        $v#6630[0] += $tmp#6629[0]
        $v#6630[1] += $tmp#6629[1]
        $v#6630[2] += $tmp#6629[2]
        $rhoVelocityFactorSkew#6626 = $v#6630
        $rhoEnergyFactorSkew#6627 += ($Fluid#6142[$7720].rhoEnthalpy*$tmp#6628)
        $tmp#6628 = (double(0.5)*$Fluid#6142[$7720].velocity[int32(1)])
        $rhoFactorSkew#6625 += ($Fluid#6142[$7721].rho*$tmp#6628)
        var $tmp#6631 : double[3] = vs_mul_double_3($Fluid#6142[$7721].rhoVelocity, $tmp#6628)
        var $v#6632 : double[3] = $rhoVelocityFactorSkew#6626
        $v#6632[0] += $tmp#6631[0]
        $v#6632[1] += $tmp#6631[1]
        $v#6632[2] += $tmp#6631[2]
        $rhoVelocityFactorSkew#6626 = $v#6632
        $rhoEnergyFactorSkew#6627 += ($Fluid#6142[$7721].rhoEnthalpy*$tmp#6628)
        var $s#6633 : double = double(0.5)
        var $rhoFlux_temp#6634 : double = (($s#6633*$rhoFactorDiagonal#6621)+((double(int32(1))-$s#6633)*$rhoFactorSkew#6625))
        var $rhoVelocityFlux_temp#6635 : double[3] = vv_add_double_3(vs_mul_double_3($rhoVelocityFactorDiagonal#6622, $s#6633), vs_mul_double_3($rhoVelocityFactorSkew#6626, (double(int32(1))-$s#6633)))
        var $rhoEnergyFlux_temp#6636 : double = (($s#6633*$rhoEnergyFactorDiagonal#6623)+((double(int32(1))-$s#6633)*$rhoEnergyFactorSkew#6627))
        $rhoVelocityFlux_temp#6635[int32(1)] += $fpdiag#6624
        $7722 = array($rhoFlux_temp#6634, $rhoVelocityFlux_temp#6635[int32(0)], $rhoVelocityFlux_temp#6635[int32(1)], $rhoVelocityFlux_temp#6635[int32(2)], $rhoEnergyFlux_temp#6636)
      end
      var $dummy#7723 : int32 = 0
      var $flux#7547 : double[5] = $7722
      $Fluid#6142[$c#7529].rhoFluxY = $flux#7547[int32(0)]
      $Fluid#6142[$c#7529].rhoVelocityFluxY = array($flux#7547[int32(1)], $flux#7547[int32(2)], $flux#7547[int32(3)])
      $Fluid#6142[$c#7529].rhoEnergyFluxY = $flux#7547[int32(4)]
      var $7725 : double
      do
        var $7724 : double = $Fluid#6142[(($c#7529+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        var $viscosity#5391 : double = double(double(0))
        if (int32(3003)==int32(3001)) then
          $viscosity#5391 = double(0.2454194)
        else
          if (int32(3003)==int32(3002)) then
            $viscosity#5391 = (double(0.9683056)*pow(($7724/double(0.9913615)), double(0.75)))
          else
            if (int32(3003)==int32(3003)) then
              $viscosity#5391 = ((double(0.5815515)*pow(($7724/double(0.6935851)), (double(3)/double(2))))*((double(0.6935851)+double(0.2772065))/($7724+double(0.2772065))))
            else
              std.assert(false, "(Liszt assertion)")
            end
          end
        end
        $7725 = $viscosity#5391
      end
      var $dummy#7726 : int32 = 0
      var $7728 : double
      do
        var $7727 : double = $Fluid#6142[$c#7529].temperature
        var $viscosity#5391 : double = double(double(0))
        if (int32(3003)==int32(3001)) then
          $viscosity#5391 = double(0.2454194)
        else
          if (int32(3003)==int32(3002)) then
            $viscosity#5391 = (double(0.9683056)*pow(($7727/double(0.9913615)), double(0.75)))
          else
            if (int32(3003)==int32(3003)) then
              $viscosity#5391 = ((double(0.5815515)*pow(($7727/double(0.6935851)), (double(3)/double(2))))*((double(0.6935851)+double(0.2772065))/($7727+double(0.2772065))))
            else
              std.assert(false, "(Liszt assertion)")
            end
          end
        end
        $7728 = $viscosity#5391
      end
      var $dummy#7729 : int32 = 0
      var $muFace#7548 : double = (double(0.5)*($7728+$7725))
      var $velocityFace#7549 : double[3] = double[3](array(double(0), double(0), double(0)))
      var $velocityY_XFace#7550 : double = double(double(0))
      var $velocityY_ZFace#7551 : double = double(double(0))
      var $velocityX_XFace#7552 : double = double(double(0))
      var $velocityZ_ZFace#7553 : double = double(double(0))
      $velocityFace#7549 = vs_mul_double_3(vv_add_double_3($Fluid#6142[$c#7529].velocity, $Fluid#6142[(($c#7529+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity), double(0.5))
      $velocityY_XFace#7550 = (double(0.5)*($Fluid#6142[$c#7529].velocityGradientX[int32(1)]+$Fluid#6142[(($c#7529+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocityGradientX[int32(1)]))
      $velocityY_ZFace#7551 = (double(0.5)*($Fluid#6142[$c#7529].velocityGradientZ[int32(1)]+$Fluid#6142[(($c#7529+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocityGradientZ[int32(1)]))
      $velocityX_XFace#7552 = (double(0.5)*($Fluid#6142[$c#7529].velocityGradientX[int32(0)]+$Fluid#6142[(($c#7529+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocityGradientX[int32(0)]))
      $velocityZ_ZFace#7553 = (double(0.5)*($Fluid#6142[$c#7529].velocityGradientZ[int32(2)]+$Fluid#6142[(($c#7529+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocityGradientZ[int32(2)]))
      var $velocityX_YFace#7554 : double = double(double(0))
      var $velocityY_YFace#7555 : double = double(double(0))
      var $velocityZ_YFace#7556 : double = double(double(0))
      var $temperature_YFace#7557 : double = double(double(0))
      $velocityX_YFace#7554 = (double(0.5)*($Fluid#6142[(($c#7529+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity[int32(0)]-$Fluid#6142[$c#7529].velocity[int32(0)]))
      $velocityY_YFace#7555 = (double(0.5)*($Fluid#6142[(($c#7529+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity[int32(1)]-$Fluid#6142[$c#7529].velocity[int32(1)]))
      $velocityZ_YFace#7556 = (double(0.5)*($Fluid#6142[(($c#7529+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity[int32(2)]-$Fluid#6142[$c#7529].velocity[int32(2)]))
      $temperature_YFace#7557 = (double(0.5)*($Fluid#6142[(($c#7529+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature-$Fluid#6142[$c#7529].temperature))
      $velocityX_YFace#7554 *= (1/(double(0.0010188725631769)*double(0.5)))
      $velocityY_YFace#7555 *= (1/(double(0.0010188725631769)*double(0.5)))
      $velocityZ_YFace#7556 *= (1/(double(0.0010188725631769)*double(0.5)))
      $temperature_YFace#7557 *= (1/(double(0.0010188725631769)*double(0.5)))
      var $sigmaXY#7558 : double = ($muFace#7548*($velocityX_YFace#7554+$velocityY_XFace#7550))
      var $sigmaYY#7559 : double = (($muFace#7548*(((double(4)*$velocityY_YFace#7555)-(double(2)*$velocityX_XFace#7552))-(double(2)*$velocityZ_ZFace#7553)))/double(3))
      var $sigmaZY#7560 : double = ($muFace#7548*($velocityZ_YFace#7556+$velocityY_ZFace#7551))
      var $usigma#7561 : double = ((($velocityFace#7549[int32(0)]*$sigmaXY#7558)+($velocityFace#7549[int32(1)]*$sigmaYY#7559))+($velocityFace#7549[int32(2)]*$sigmaZY#7560))
      var $cp#7562 : double = ((double(0.8198141)*double(0.7272462))/(double(0.8198141)-double(1)))
      var $heatFlux#7563 : double = ((-(($cp#7562*$muFace#7548)/double(0.1565651)))*$temperature_YFace#7557)
      $Fluid#6142[$c#7529].rhoVelocityFluxY[int32(0)] += (-$sigmaXY#7558)
      $Fluid#6142[$c#7529].rhoVelocityFluxY[int32(1)] += (-$sigmaYY#7559)
      $Fluid#6142[$c#7529].rhoVelocityFluxY[int32(2)] += (-$sigmaZY#7560)
      $Fluid#6142[$c#7529].rhoEnergyFluxY += (-($usigma#7561-$heatFlux#7563))
    else
    end
    if ((not ((((((max(int32((uint64(int32(0))-int3d($c#7529).x)), int32(0))>int32(0)) or (max(int32((int3d($c#7529).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#7529).y)), int32(0))>int32(0))) or (max(int32((int3d($c#7529).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#7529).z)), int32(0))>int32(0))) or (max(int32((int3d($c#7529).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0)))) or (max(int32((uint64(int32(1))-int3d($c#7529).z)), int32(0))==int32(1))) then
      var $7732 : double[5]
      do
        var $7730 : int3d = int3d($c#7529)
        var $7731 : int3d = (($c#7529+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))
        var $rhoFactorDiagonal#7041 : double = double(double(0))
        var $rhoVelocityFactorDiagonal#7042 : double[3] = double[3](array(double(0), double(0), double(0)))
        var $rhoEnergyFactorDiagonal#7043 : double = double(double(0))
        var $fpdiag#7044 : double = double(double(0))
        $rhoFactorDiagonal#7041 = (double(0.5)*(($Fluid#6142[$7730].rho*$Fluid#6142[$7730].velocity[int32(2)])+($Fluid#6142[$7731].rho*$Fluid#6142[$7731].velocity[int32(2)])))
        $rhoVelocityFactorDiagonal#7042 = vs_mul_double_3(vv_add_double_3(vs_mul_double_3($Fluid#6142[$7730].rhoVelocity, $Fluid#6142[$7730].velocity[int32(2)]), vs_mul_double_3($Fluid#6142[$7731].rhoVelocity, $Fluid#6142[$7731].velocity[int32(2)])), double(0.5))
        $rhoEnergyFactorDiagonal#7043 = (double(0.5)*(($Fluid#6142[$7730].rhoEnthalpy*$Fluid#6142[$7730].velocity[int32(2)])+($Fluid#6142[$7731].rhoEnthalpy*$Fluid#6142[$7731].velocity[int32(2)])))
        $fpdiag#7044 += (double(0.5)*($Fluid#6142[$7730].pressure+$Fluid#6142[$7731].pressure))
        var $rhoFactorSkew#7045 : double = double(double(0))
        var $rhoVelocityFactorSkew#7046 : double[3] = double[3](array(double(0), double(0), double(0)))
        var $rhoEnergyFactorSkew#7047 : double = double(double(0))
        var $tmp#7048 : double = double(double(0))
        $tmp#7048 = (double(0.5)*$Fluid#6142[$7731].velocity[int32(2)])
        $rhoFactorSkew#7045 += ($Fluid#6142[$7730].rho*$tmp#7048)
        var $tmp#7049 : double[3] = vs_mul_double_3($Fluid#6142[$7730].rhoVelocity, $tmp#7048)
        var $v#7050 : double[3] = $rhoVelocityFactorSkew#7046
        $v#7050[0] += $tmp#7049[0]
        $v#7050[1] += $tmp#7049[1]
        $v#7050[2] += $tmp#7049[2]
        $rhoVelocityFactorSkew#7046 = $v#7050
        $rhoEnergyFactorSkew#7047 += ($Fluid#6142[$7730].rhoEnthalpy*$tmp#7048)
        $tmp#7048 = (double(0.5)*$Fluid#6142[$7730].velocity[int32(2)])
        $rhoFactorSkew#7045 += ($Fluid#6142[$7731].rho*$tmp#7048)
        var $tmp#7051 : double[3] = vs_mul_double_3($Fluid#6142[$7731].rhoVelocity, $tmp#7048)
        var $v#7052 : double[3] = $rhoVelocityFactorSkew#7046
        $v#7052[0] += $tmp#7051[0]
        $v#7052[1] += $tmp#7051[1]
        $v#7052[2] += $tmp#7051[2]
        $rhoVelocityFactorSkew#7046 = $v#7052
        $rhoEnergyFactorSkew#7047 += ($Fluid#6142[$7731].rhoEnthalpy*$tmp#7048)
        var $s#7053 : double = double(0.5)
        var $rhoFlux_temp#7054 : double = (($s#7053*$rhoFactorDiagonal#7041)+((double(int32(1))-$s#7053)*$rhoFactorSkew#7045))
        var $rhoVelocityFlux_temp#7055 : double[3] = vv_add_double_3(vs_mul_double_3($rhoVelocityFactorDiagonal#7042, $s#7053), vs_mul_double_3($rhoVelocityFactorSkew#7046, (double(int32(1))-$s#7053)))
        var $rhoEnergyFlux_temp#7056 : double = (($s#7053*$rhoEnergyFactorDiagonal#7043)+((double(int32(1))-$s#7053)*$rhoEnergyFactorSkew#7047))
        $rhoVelocityFlux_temp#7055[int32(2)] += $fpdiag#7044
        $7732 = array($rhoFlux_temp#7054, $rhoVelocityFlux_temp#7055[int32(0)], $rhoVelocityFlux_temp#7055[int32(1)], $rhoVelocityFlux_temp#7055[int32(2)], $rhoEnergyFlux_temp#7056)
      end
      var $dummy#7733 : int32 = 0
      var $flux#7564 : double[5] = $7732
      $Fluid#6142[$c#7529].rhoFluxZ = $flux#7564[int32(0)]
      $Fluid#6142[$c#7529].rhoVelocityFluxZ = array($flux#7564[int32(1)], $flux#7564[int32(2)], $flux#7564[int32(3)])
      $Fluid#6142[$c#7529].rhoEnergyFluxZ = $flux#7564[int32(4)]
      var $7735 : double
      do
        var $7734 : double = $Fluid#6142[(($c#7529+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        var $viscosity#5391 : double = double(double(0))
        if (int32(3003)==int32(3001)) then
          $viscosity#5391 = double(0.2454194)
        else
          if (int32(3003)==int32(3002)) then
            $viscosity#5391 = (double(0.9683056)*pow(($7734/double(0.9913615)), double(0.75)))
          else
            if (int32(3003)==int32(3003)) then
              $viscosity#5391 = ((double(0.5815515)*pow(($7734/double(0.6935851)), (double(3)/double(2))))*((double(0.6935851)+double(0.2772065))/($7734+double(0.2772065))))
            else
              std.assert(false, "(Liszt assertion)")
            end
          end
        end
        $7735 = $viscosity#5391
      end
      var $dummy#7736 : int32 = 0
      var $7738 : double
      do
        var $7737 : double = $Fluid#6142[$c#7529].temperature
        var $viscosity#5391 : double = double(double(0))
        if (int32(3003)==int32(3001)) then
          $viscosity#5391 = double(0.2454194)
        else
          if (int32(3003)==int32(3002)) then
            $viscosity#5391 = (double(0.9683056)*pow(($7737/double(0.9913615)), double(0.75)))
          else
            if (int32(3003)==int32(3003)) then
              $viscosity#5391 = ((double(0.5815515)*pow(($7737/double(0.6935851)), (double(3)/double(2))))*((double(0.6935851)+double(0.2772065))/($7737+double(0.2772065))))
            else
              std.assert(false, "(Liszt assertion)")
            end
          end
        end
        $7738 = $viscosity#5391
      end
      var $dummy#7739 : int32 = 0
      var $muFace#7565 : double = (double(0.5)*($7738+$7735))
      var $velocityFace#7566 : double[3] = double[3](array(double(0), double(0), double(0)))
      var $velocityZ_XFace#7567 : double = double(double(0))
      var $velocityZ_YFace#7568 : double = double(double(0))
      var $velocityX_XFace#7569 : double = double(double(0))
      var $velocityY_YFace#7570 : double = double(double(0))
      $velocityFace#7566 = vs_mul_double_3(vv_add_double_3($Fluid#6142[$c#7529].velocity, $Fluid#6142[(($c#7529+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity), double(0.5))
      $velocityZ_XFace#7567 = (double(0.5)*($Fluid#6142[$c#7529].velocityGradientX[int32(2)]+$Fluid#6142[(($c#7529+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocityGradientX[int32(2)]))
      $velocityZ_YFace#7568 = (double(0.5)*($Fluid#6142[$c#7529].velocityGradientY[int32(2)]+$Fluid#6142[(($c#7529+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocityGradientY[int32(2)]))
      $velocityX_XFace#7569 = (double(0.5)*($Fluid#6142[$c#7529].velocityGradientX[int32(0)]+$Fluid#6142[(($c#7529+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocityGradientX[int32(0)]))
      $velocityY_YFace#7570 = (double(0.5)*($Fluid#6142[$c#7529].velocityGradientY[int32(1)]+$Fluid#6142[(($c#7529+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocityGradientY[int32(1)]))
      var $velocityX_ZFace#7571 : double = double(double(0))
      var $velocityY_ZFace#7572 : double = double(double(0))
      var $velocityZ_ZFace#7573 : double = double(double(0))
      var $temperature_ZFace#7574 : double = double(double(0))
      $velocityX_ZFace#7571 = (double(0.5)*($Fluid#6142[(($c#7529+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity[int32(0)]-$Fluid#6142[$c#7529].velocity[int32(0)]))
      $velocityY_ZFace#7572 = (double(0.5)*($Fluid#6142[(($c#7529+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity[int32(1)]-$Fluid#6142[$c#7529].velocity[int32(1)]))
      $velocityZ_ZFace#7573 = (double(0.5)*($Fluid#6142[(($c#7529+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity[int32(2)]-$Fluid#6142[$c#7529].velocity[int32(2)]))
      $temperature_ZFace#7574 = (double(0.5)*($Fluid#6142[(($c#7529+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature-$Fluid#6142[$c#7529].temperature))
      $velocityX_ZFace#7571 *= (1/(double(0.0011352949453552)*double(0.5)))
      $velocityY_ZFace#7572 *= (1/(double(0.0011352949453552)*double(0.5)))
      $velocityZ_ZFace#7573 *= (1/(double(0.0011352949453552)*double(0.5)))
      $temperature_ZFace#7574 *= (1/(double(0.0011352949453552)*double(0.5)))
      var $sigmaXZ#7575 : double = ($muFace#7565*($velocityX_ZFace#7571+$velocityZ_XFace#7567))
      var $sigmaYZ#7576 : double = ($muFace#7565*($velocityY_ZFace#7572+$velocityZ_YFace#7568))
      var $sigmaZZ#7577 : double = (($muFace#7565*(((double(4)*$velocityZ_ZFace#7573)-(double(2)*$velocityX_XFace#7569))-(double(2)*$velocityY_YFace#7570)))/double(3))
      var $usigma#7578 : double = ((($velocityFace#7566[int32(0)]*$sigmaXZ#7575)+($velocityFace#7566[int32(1)]*$sigmaYZ#7576))+($velocityFace#7566[int32(2)]*$sigmaZZ#7577))
      var $cp#7579 : double = ((double(0.8198141)*double(0.7272462))/(double(0.8198141)-double(1)))
      var $heatFlux#7580 : double = ((-(($cp#7579*$muFace#7565)/double(0.1565651)))*$temperature_ZFace#7574)
      $Fluid#6142[$c#7529].rhoVelocityFluxZ[int32(0)] += (-$sigmaXZ#7575)
      $Fluid#6142[$c#7529].rhoVelocityFluxZ[int32(1)] += (-$sigmaYZ#7576)
      $Fluid#6142[$c#7529].rhoVelocityFluxZ[int32(2)] += (-$sigmaZZ#7577)
      $Fluid#6142[$c#7529].rhoEnergyFluxZ += (-($usigma#7578-$heatFlux#7580))
    else
    end
    var $v#7581 : int32 = int32(0)
    if ($v#7581==int32(1)) then
      var $tmp1#7582 : double = $Fluid#6142[(($c#7529+{-1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity[int32(0)]
      var $tmp2#7583 : double = $Fluid#6142[(($c#7529+{0, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity[int32(1)]
      var $tmp3#7584 : double = $Fluid#6142[(($c#7529+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity[int32(2)]
    else
    end
  end
end
task Flow_AddUpdateUsingFlux($dom#8390 : region#2121(ispace#2121(int3d), Fluid_columns), $Fluid#8392 : region#2122(ispace#2122(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#8392.rhoEnergyFluxX), reads($Fluid#8392.rhoEnergyFluxY), reads($Fluid#8392.rhoEnergyFluxZ), reads($Fluid#8392.rhoEnergy_t), writes($Fluid#8392.rhoEnergy_t), reads($Fluid#8392.rhoFluxX), reads($Fluid#8392.rhoFluxY), reads($Fluid#8392.rhoFluxZ), reads($Fluid#8392.rhoVelocityFluxX), reads($Fluid#8392.rhoVelocityFluxY), reads($Fluid#8392.rhoVelocityFluxZ), reads($Fluid#8392.rhoVelocity_t), writes($Fluid#8392.rhoVelocity_t), reads($Fluid#8392.rho_t), writes($Fluid#8392.rho_t), $dom#8390 <= $Fluid#8392
do
  for $c#8437 : int3d(Fluid_columns, $dom#8390) in $dom#8390 do
    if (not ((((((max(int32((uint64(int32(0))-int3d($c#8437).x)), int32(0))>int32(0)) or (max(int32((int3d($c#8437).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#8437).y)), int32(0))>int32(0))) or (max(int32((int3d($c#8437).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#8437).z)), int32(0))>int32(0))) or (max(int32((int3d($c#8437).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0)))) then
      $Fluid#8392[$c#8437].rho_t += ((-($Fluid#8392[$c#8437].rhoFluxX-$Fluid#8392[(($c#8437+{-1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].rhoFluxX))/double(0.0015841600529101))
      var $tmp#8438 : double[3] = vs_div_double_3(vs_mul_double_3(vv_sub_double_3($Fluid#8392[$c#8437].rhoVelocityFluxX, $Fluid#8392[(($c#8437+{-1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].rhoVelocityFluxX), double((-1))), double(0.0015841600529101))
      var $v#8439 : double[3] = $Fluid#8392[$c#8437].rhoVelocity_t
      $v#8439[0] += $tmp#8438[0]
      $v#8439[1] += $tmp#8438[1]
      $v#8439[2] += $tmp#8438[2]
      $Fluid#8392[$c#8437].rhoVelocity_t = $v#8439
      $Fluid#8392[$c#8437].rhoEnergy_t += ((-($Fluid#8392[$c#8437].rhoEnergyFluxX-$Fluid#8392[(($c#8437+{-1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].rhoEnergyFluxX))/double(0.0015841600529101))
      $Fluid#8392[$c#8437].rho_t += ((-($Fluid#8392[$c#8437].rhoFluxY-$Fluid#8392[(($c#8437+{0, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].rhoFluxY))/double(0.0010188725631769))
      var $tmp#8440 : double[3] = vs_div_double_3(vs_mul_double_3(vv_sub_double_3($Fluid#8392[$c#8437].rhoVelocityFluxY, $Fluid#8392[(($c#8437+{0, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].rhoVelocityFluxY), double((-1))), double(0.0010188725631769))
      var $v#8441 : double[3] = $Fluid#8392[$c#8437].rhoVelocity_t
      $v#8441[0] += $tmp#8440[0]
      $v#8441[1] += $tmp#8440[1]
      $v#8441[2] += $tmp#8440[2]
      $Fluid#8392[$c#8437].rhoVelocity_t = $v#8441
      $Fluid#8392[$c#8437].rhoEnergy_t += ((-($Fluid#8392[$c#8437].rhoEnergyFluxY-$Fluid#8392[(($c#8437+{0, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].rhoEnergyFluxY))/double(0.0010188725631769))
      $Fluid#8392[$c#8437].rho_t += ((-($Fluid#8392[$c#8437].rhoFluxZ-$Fluid#8392[(($c#8437+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].rhoFluxZ))/double(0.0011352949453552))
      var $tmp#8442 : double[3] = vs_div_double_3(vs_mul_double_3(vv_sub_double_3($Fluid#8392[$c#8437].rhoVelocityFluxZ, $Fluid#8392[(($c#8437+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].rhoVelocityFluxZ), double((-1))), double(0.0011352949453552))
      var $v#8443 : double[3] = $Fluid#8392[$c#8437].rhoVelocity_t
      $v#8443[0] += $tmp#8442[0]
      $v#8443[1] += $tmp#8442[1]
      $v#8443[2] += $tmp#8442[2]
      $Fluid#8392[$c#8437].rhoVelocity_t = $v#8443
      $Fluid#8392[$c#8437].rhoEnergy_t += ((-($Fluid#8392[$c#8437].rhoEnergyFluxZ-$Fluid#8392[(($c#8437+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].rhoEnergyFluxZ))/double(0.0011352949453552))
    else
    end
  end
end
task Flow_AddBodyForces($dom#8590 : region#2177(ispace#2177(int3d), Fluid_columns), $Fluid#8592 : region#2178(ispace#2178(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#8592.rho), reads($Fluid#8592.rhoEnergy_t), writes($Fluid#8592.rhoEnergy_t), reads($Fluid#8592.rhoVelocity_t), writes($Fluid#8592.rhoVelocity_t), reads($Fluid#8592.velocity), $dom#8590 <= $Fluid#8592
do
  for $c#8609 : int3d(Fluid_columns, $dom#8590) in $dom#8590 do
    if (not ((((((max(int32((uint64(int32(0))-int3d($c#8609).x)), int32(0))>int32(0)) or (max(int32((int3d($c#8609).x-uint64(((int32(378)-int32(1))-int32(0))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#8609).y)), int32(0))>int32(0))) or (max(int32((int3d($c#8609).y-uint64(((int32(833)-int32(1))-int32(1))))), int32(0))>int32(0))) or (max(int32((uint64(int32(1))-int3d($c#8609).z)), int32(0))>int32(0))) or (max(int32((int3d($c#8609).z-uint64(((int32(734)-int32(1))-int32(1))))), int32(0))>int32(0)))) then
      var $tmp#8610 : double[3] = vs_mul_double_3(array(double(0.7751856), double(0.7816449), double(0.1264908)), $Fluid#8592[$c#8609].rho)
      var $v#8611 : double[3] = $Fluid#8592[$c#8609].rhoVelocity_t
      $v#8611[0] += $tmp#8610[0]
      $v#8611[1] += $tmp#8610[1]
      $v#8611[2] += $tmp#8610[2]
      $Fluid#8592[$c#8609].rhoVelocity_t = $v#8611
      $Fluid#8592[$c#8609].rhoEnergy_t += ($Fluid#8592[$c#8609].rho*dot_double_3(array(double(0.7751856), double(0.7816449), double(0.1264908)), $Fluid#8592[$c#8609].velocity))
    else
    end
  end
end
task Particles_LocateInCells($dom#8661 : region#2192(ispace#2192(int1d), particles_columns), $particles#8663 : region#2193(ispace#2193(int1d), particles_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($particles#8663.cell), writes($particles#8663.cell), reads($particles#8663.position), reads($particles#8663.__valid), $dom#8661 <= $particles#8663
do
  for $p#8666 : int1d(particles_columns, $dom#8661) in $dom#8661 do
    if $particles#8663[$p#8666].__valid then
      $particles#8663[$p#8666].cell = int3d({(uint64((fmod(((($particles#8663[$p#8666].position[int32(0)]-double(0.6987052))/double(0.5988125))*double(int32(378))), double(int32(378)))+double(int32(378))))%uint64(int32(378))), uint64(max(double(0), min(double((int32(833)-int32(1))), ((($particles#8663[$p#8666].position[int32(1)]-double(0.99700342743682))/double(0.84872084512635))*double(int32(833)))))), uint64(max(double(0), min(double((int32(734)-int32(1))), ((($particles#8663[$p#8666].position[int32(2)]-double(0.70456500505464))/double(0.83330648989071))*double(int32(734))))))})
    else
    end
  end
end
task Fluid_elemColor($idx#8693 : int3d) : int3d
-- leaf (false), inner (false), idempotent (false)
  $idx#8693.x = min(max($idx#8693.x, 0), ((378+0)-1))
  $idx#8693.y = min(max($idx#8693.y, 1), ((833+1)-1))
  $idx#8693.z = min(max($idx#8693.z, 1), ((734+1)-1))
  return int3d({(($idx#8693.x-0)/189), (($idx#8693.y-1)/277), (($idx#8693.z-1)/183)})
end
terra particles_pushElement(dst : &opaque,idx : int32,src : particles_columns) : {}
    var ptr : &int8 = [&int8](dst) + idx * 376
    memcpy([&opaque](ptr), [&opaque](&src), [uint64](376))
end
terra particles_getBasePointer(pr : legion_physical_region_t,fid : uint32,runtime : legion_runtime_t) : &opaque
    var acc : legion_accessor_array_1d_t = legion_physical_region_get_field_accessor_array_1d(pr, fid)
    var lr : legion_logical_region_t = legion_physical_region_get_logical_region(pr)
    var domain : legion_domain_t = legion_index_space_get_domain(runtime, lr.index_space)
    var rect : legion_rect_1d_t = legion_domain_get_rect_1d(domain)
    var subrect : legion_rect_1d_t
    var offsets : legion_byte_offset_t[1]
    var p : &opaque = legion_accessor_array_1d_raw_rect_ptr(acc, rect, &subrect, &[&legion_byte_offset_t](offsets)[0])
    legion_accessor_array_1d_destroy(acc)
    return p
end
terra particles_getOffset() : int64
    var x : particles_columns
    return [&int8](&x.__valid) - [&int8](&x)
end
task particles_pushAll($partColor#8698 : int3d, $r#8696 : region#2204(ispace#2204(int1d), particles_columns), $q0#8702 : region#2206(ispace#2206(int1d), int8[376]), $q1#8714 : region#2207(ispace#2207(int1d), int8[376]), $q2#8726 : region#2208(ispace#2208(int1d), int8[376]), $q3#8738 : region#2209(ispace#2209(int1d), int8[376]), $q4#8750 : region#2210(ispace#2210(int1d), int8[376]), $q5#8762 : region#2211(ispace#2211(int1d), int8[376]), $q6#8774 : region#2212(ispace#2212(int1d), int8[376]), $q7#8786 : region#2213(ispace#2213(int1d), int8[376]), $q8#8798 : region#2214(ispace#2214(int1d), int8[376]), $q9#8810 : region#2215(ispace#2215(int1d), int8[376]), $q10#8822 : region#2216(ispace#2216(int1d), int8[376]), $q11#8834 : region#2217(ispace#2217(int1d), int8[376]), $q12#8846 : region#2218(ispace#2218(int1d), int8[376]), $q13#8858 : region#2219(ispace#2219(int1d), int8[376]), $q14#8870 : region#2220(ispace#2220(int1d), int8[376]), $q15#8882 : region#2221(ispace#2221(int1d), int8[376]), $q16#8894 : region#2222(ispace#2222(int1d), int8[376]), $q17#8906 : region#2223(ispace#2223(int1d), int8[376]), $q18#8918 : region#2224(ispace#2224(int1d), int8[376]), $q19#8930 : region#2225(ispace#2225(int1d), int8[376]), $q20#8942 : region#2226(ispace#2226(int1d), int8[376]), $q21#8954 : region#2227(ispace#2227(int1d), int8[376]), $q22#8966 : region#2228(ispace#2228(int1d), int8[376]), $q23#8978 : region#2229(ispace#2229(int1d), int8[376]), $q24#8990 : region#2230(ispace#2230(int1d), int8[376]), $q25#9002 : region#2231(ispace#2231(int1d), int8[376]))
-- leaf (false), inner (false), idempotent (false)
where
  reads($r#8696), writes($r#8696.__valid), reads($q0#8702), writes($q0#8702), reads($q1#8714), writes($q1#8714), reads($q2#8726), writes($q2#8726), reads($q3#8738), writes($q3#8738), reads($q4#8750), writes($q4#8750), reads($q5#8762), writes($q5#8762), reads($q6#8774), writes($q6#8774), reads($q7#8786), writes($q7#8786), reads($q8#8798), writes($q8#8798), reads($q9#8810), writes($q9#8810), reads($q10#8822), writes($q10#8822), reads($q11#8834), writes($q11#8834), reads($q12#8846), writes($q12#8846), reads($q13#8858), writes($q13#8858), reads($q14#8870), writes($q14#8870), reads($q15#8882), writes($q15#8882), reads($q16#8894), writes($q16#8894), reads($q17#8906), writes($q17#8906), reads($q18#8918), writes($q18#8918), reads($q19#8930), writes($q19#8930), reads($q20#8942), writes($q20#8942), reads($q21#8954), writes($q21#8954), reads($q22#8966), writes($q22#8966), reads($q23#8978), writes($q23#8978), reads($q24#8990), writes($q24#8990), reads($q25#9002), writes($q25#9002)
do
  for $qPtr#9013 : int1d(int8[376], $q0#8702) in $q0#8702 do
    $q0#8702[$qPtr#9013][368LL] = int8(false)
  end
  var $qBasePtr0#9014 : &opaque = particles_getBasePointer(__physical($q0#8702)[0], __fields($q0#8702)[0], __runtime())
  for $qPtr#9015 : int1d(int8[376], $q1#8714) in $q1#8714 do
    $q1#8714[$qPtr#9015][368LL] = int8(false)
  end
  var $qBasePtr1#9016 : &opaque = particles_getBasePointer(__physical($q1#8714)[0], __fields($q1#8714)[0], __runtime())
  for $qPtr#9017 : int1d(int8[376], $q2#8726) in $q2#8726 do
    $q2#8726[$qPtr#9017][368LL] = int8(false)
  end
  var $qBasePtr2#9018 : &opaque = particles_getBasePointer(__physical($q2#8726)[0], __fields($q2#8726)[0], __runtime())
  for $qPtr#9019 : int1d(int8[376], $q3#8738) in $q3#8738 do
    $q3#8738[$qPtr#9019][368LL] = int8(false)
  end
  var $qBasePtr3#9020 : &opaque = particles_getBasePointer(__physical($q3#8738)[0], __fields($q3#8738)[0], __runtime())
  for $qPtr#9021 : int1d(int8[376], $q4#8750) in $q4#8750 do
    $q4#8750[$qPtr#9021][368LL] = int8(false)
  end
  var $qBasePtr4#9022 : &opaque = particles_getBasePointer(__physical($q4#8750)[0], __fields($q4#8750)[0], __runtime())
  for $qPtr#9023 : int1d(int8[376], $q5#8762) in $q5#8762 do
    $q5#8762[$qPtr#9023][368LL] = int8(false)
  end
  var $qBasePtr5#9024 : &opaque = particles_getBasePointer(__physical($q5#8762)[0], __fields($q5#8762)[0], __runtime())
  for $qPtr#9025 : int1d(int8[376], $q6#8774) in $q6#8774 do
    $q6#8774[$qPtr#9025][368LL] = int8(false)
  end
  var $qBasePtr6#9026 : &opaque = particles_getBasePointer(__physical($q6#8774)[0], __fields($q6#8774)[0], __runtime())
  for $qPtr#9027 : int1d(int8[376], $q7#8786) in $q7#8786 do
    $q7#8786[$qPtr#9027][368LL] = int8(false)
  end
  var $qBasePtr7#9028 : &opaque = particles_getBasePointer(__physical($q7#8786)[0], __fields($q7#8786)[0], __runtime())
  for $qPtr#9029 : int1d(int8[376], $q8#8798) in $q8#8798 do
    $q8#8798[$qPtr#9029][368LL] = int8(false)
  end
  var $qBasePtr8#9030 : &opaque = particles_getBasePointer(__physical($q8#8798)[0], __fields($q8#8798)[0], __runtime())
  for $qPtr#9031 : int1d(int8[376], $q9#8810) in $q9#8810 do
    $q9#8810[$qPtr#9031][368LL] = int8(false)
  end
  var $qBasePtr9#9032 : &opaque = particles_getBasePointer(__physical($q9#8810)[0], __fields($q9#8810)[0], __runtime())
  for $qPtr#9033 : int1d(int8[376], $q10#8822) in $q10#8822 do
    $q10#8822[$qPtr#9033][368LL] = int8(false)
  end
  var $qBasePtr10#9034 : &opaque = particles_getBasePointer(__physical($q10#8822)[0], __fields($q10#8822)[0], __runtime())
  for $qPtr#9035 : int1d(int8[376], $q11#8834) in $q11#8834 do
    $q11#8834[$qPtr#9035][368LL] = int8(false)
  end
  var $qBasePtr11#9036 : &opaque = particles_getBasePointer(__physical($q11#8834)[0], __fields($q11#8834)[0], __runtime())
  for $qPtr#9037 : int1d(int8[376], $q12#8846) in $q12#8846 do
    $q12#8846[$qPtr#9037][368LL] = int8(false)
  end
  var $qBasePtr12#9038 : &opaque = particles_getBasePointer(__physical($q12#8846)[0], __fields($q12#8846)[0], __runtime())
  for $qPtr#9039 : int1d(int8[376], $q13#8858) in $q13#8858 do
    $q13#8858[$qPtr#9039][368LL] = int8(false)
  end
  var $qBasePtr13#9040 : &opaque = particles_getBasePointer(__physical($q13#8858)[0], __fields($q13#8858)[0], __runtime())
  for $qPtr#9041 : int1d(int8[376], $q14#8870) in $q14#8870 do
    $q14#8870[$qPtr#9041][368LL] = int8(false)
  end
  var $qBasePtr14#9042 : &opaque = particles_getBasePointer(__physical($q14#8870)[0], __fields($q14#8870)[0], __runtime())
  for $qPtr#9043 : int1d(int8[376], $q15#8882) in $q15#8882 do
    $q15#8882[$qPtr#9043][368LL] = int8(false)
  end
  var $qBasePtr15#9044 : &opaque = particles_getBasePointer(__physical($q15#8882)[0], __fields($q15#8882)[0], __runtime())
  for $qPtr#9045 : int1d(int8[376], $q16#8894) in $q16#8894 do
    $q16#8894[$qPtr#9045][368LL] = int8(false)
  end
  var $qBasePtr16#9046 : &opaque = particles_getBasePointer(__physical($q16#8894)[0], __fields($q16#8894)[0], __runtime())
  for $qPtr#9047 : int1d(int8[376], $q17#8906) in $q17#8906 do
    $q17#8906[$qPtr#9047][368LL] = int8(false)
  end
  var $qBasePtr17#9048 : &opaque = particles_getBasePointer(__physical($q17#8906)[0], __fields($q17#8906)[0], __runtime())
  for $qPtr#9049 : int1d(int8[376], $q18#8918) in $q18#8918 do
    $q18#8918[$qPtr#9049][368LL] = int8(false)
  end
  var $qBasePtr18#9050 : &opaque = particles_getBasePointer(__physical($q18#8918)[0], __fields($q18#8918)[0], __runtime())
  for $qPtr#9051 : int1d(int8[376], $q19#8930) in $q19#8930 do
    $q19#8930[$qPtr#9051][368LL] = int8(false)
  end
  var $qBasePtr19#9052 : &opaque = particles_getBasePointer(__physical($q19#8930)[0], __fields($q19#8930)[0], __runtime())
  for $qPtr#9053 : int1d(int8[376], $q20#8942) in $q20#8942 do
    $q20#8942[$qPtr#9053][368LL] = int8(false)
  end
  var $qBasePtr20#9054 : &opaque = particles_getBasePointer(__physical($q20#8942)[0], __fields($q20#8942)[0], __runtime())
  for $qPtr#9055 : int1d(int8[376], $q21#8954) in $q21#8954 do
    $q21#8954[$qPtr#9055][368LL] = int8(false)
  end
  var $qBasePtr21#9056 : &opaque = particles_getBasePointer(__physical($q21#8954)[0], __fields($q21#8954)[0], __runtime())
  for $qPtr#9057 : int1d(int8[376], $q22#8966) in $q22#8966 do
    $q22#8966[$qPtr#9057][368LL] = int8(false)
  end
  var $qBasePtr22#9058 : &opaque = particles_getBasePointer(__physical($q22#8966)[0], __fields($q22#8966)[0], __runtime())
  for $qPtr#9059 : int1d(int8[376], $q23#8978) in $q23#8978 do
    $q23#8978[$qPtr#9059][368LL] = int8(false)
  end
  var $qBasePtr23#9060 : &opaque = particles_getBasePointer(__physical($q23#8978)[0], __fields($q23#8978)[0], __runtime())
  for $qPtr#9061 : int1d(int8[376], $q24#8990) in $q24#8990 do
    $q24#8990[$qPtr#9061][368LL] = int8(false)
  end
  var $qBasePtr24#9062 : &opaque = particles_getBasePointer(__physical($q24#8990)[0], __fields($q24#8990)[0], __runtime())
  for $qPtr#9063 : int1d(int8[376], $q25#9002) in $q25#9002 do
    $q25#9002[$qPtr#9063][368LL] = int8(false)
  end
  var $qBasePtr25#9064 : &opaque = particles_getBasePointer(__physical($q25#9002)[0], __fields($q25#9002)[0], __runtime())
  for $rPtr#9065 : int1d(particles_columns, $r#8696) in $r#8696 do
    if (@$rPtr#9065).__valid then
      var $9751 : int3d
      do
        var $9750 : int3d = (@$rPtr#9065).cell
        $9750.x = min(max($9750.x, 0), ((378+0)-1))
        $9750.y = min(max($9750.y, 1), ((833+1)-1))
        $9750.z = min(max($9750.z, 1), ((734+1)-1))
        $9751 = int3d({(($9750.x-0)/189), (($9750.y-1)/277), (($9750.z-1)/183)})
      end
      var $dummy#9752 : int32 = 0
      var $elemColor#9066 : int3d = $9751
      if ($elemColor#9066~=$partColor#8698) then
        do
          var $colorOff#9067 : int3d = int3d({0, 0, 1})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9067)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9068 : int32 = 0
            for $qPtr#9069 : int1d(int8[376], $q0#8702) in $q0#8702 do
              if (not bool($q0#8702[$qPtr#9069][368LL])) then
                particles_pushElement($qBasePtr0#9014, $idx#9068, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q0#8702[$qPtr#9069][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9068 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9070 : int3d = int3d({0, 0, -1})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9070)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9071 : int32 = 0
            for $qPtr#9072 : int1d(int8[376], $q1#8714) in $q1#8714 do
              if (not bool($q1#8714[$qPtr#9072][368LL])) then
                particles_pushElement($qBasePtr1#9016, $idx#9071, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q1#8714[$qPtr#9072][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9071 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9073 : int3d = int3d({0, 1, 0})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9073)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9074 : int32 = 0
            for $qPtr#9075 : int1d(int8[376], $q2#8726) in $q2#8726 do
              if (not bool($q2#8726[$qPtr#9075][368LL])) then
                particles_pushElement($qBasePtr2#9018, $idx#9074, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q2#8726[$qPtr#9075][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9074 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9076 : int3d = int3d({0, 1, 1})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9076)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9077 : int32 = 0
            for $qPtr#9078 : int1d(int8[376], $q3#8738) in $q3#8738 do
              if (not bool($q3#8738[$qPtr#9078][368LL])) then
                particles_pushElement($qBasePtr3#9020, $idx#9077, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q3#8738[$qPtr#9078][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9077 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9079 : int3d = int3d({0, 1, -1})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9079)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9080 : int32 = 0
            for $qPtr#9081 : int1d(int8[376], $q4#8750) in $q4#8750 do
              if (not bool($q4#8750[$qPtr#9081][368LL])) then
                particles_pushElement($qBasePtr4#9022, $idx#9080, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q4#8750[$qPtr#9081][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9080 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9082 : int3d = int3d({0, -1, 0})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9082)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9083 : int32 = 0
            for $qPtr#9084 : int1d(int8[376], $q5#8762) in $q5#8762 do
              if (not bool($q5#8762[$qPtr#9084][368LL])) then
                particles_pushElement($qBasePtr5#9024, $idx#9083, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q5#8762[$qPtr#9084][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9083 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9085 : int3d = int3d({0, -1, 1})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9085)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9086 : int32 = 0
            for $qPtr#9087 : int1d(int8[376], $q6#8774) in $q6#8774 do
              if (not bool($q6#8774[$qPtr#9087][368LL])) then
                particles_pushElement($qBasePtr6#9026, $idx#9086, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q6#8774[$qPtr#9087][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9086 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9088 : int3d = int3d({0, -1, -1})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9088)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9089 : int32 = 0
            for $qPtr#9090 : int1d(int8[376], $q7#8786) in $q7#8786 do
              if (not bool($q7#8786[$qPtr#9090][368LL])) then
                particles_pushElement($qBasePtr7#9028, $idx#9089, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q7#8786[$qPtr#9090][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9089 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9091 : int3d = int3d({1, 0, 0})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9091)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9092 : int32 = 0
            for $qPtr#9093 : int1d(int8[376], $q8#8798) in $q8#8798 do
              if (not bool($q8#8798[$qPtr#9093][368LL])) then
                particles_pushElement($qBasePtr8#9030, $idx#9092, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q8#8798[$qPtr#9093][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9092 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9094 : int3d = int3d({1, 0, 1})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9094)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9095 : int32 = 0
            for $qPtr#9096 : int1d(int8[376], $q9#8810) in $q9#8810 do
              if (not bool($q9#8810[$qPtr#9096][368LL])) then
                particles_pushElement($qBasePtr9#9032, $idx#9095, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q9#8810[$qPtr#9096][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9095 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9097 : int3d = int3d({1, 0, -1})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9097)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9098 : int32 = 0
            for $qPtr#9099 : int1d(int8[376], $q10#8822) in $q10#8822 do
              if (not bool($q10#8822[$qPtr#9099][368LL])) then
                particles_pushElement($qBasePtr10#9034, $idx#9098, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q10#8822[$qPtr#9099][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9098 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9100 : int3d = int3d({1, 1, 0})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9100)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9101 : int32 = 0
            for $qPtr#9102 : int1d(int8[376], $q11#8834) in $q11#8834 do
              if (not bool($q11#8834[$qPtr#9102][368LL])) then
                particles_pushElement($qBasePtr11#9036, $idx#9101, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q11#8834[$qPtr#9102][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9101 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9103 : int3d = int3d({1, 1, 1})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9103)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9104 : int32 = 0
            for $qPtr#9105 : int1d(int8[376], $q12#8846) in $q12#8846 do
              if (not bool($q12#8846[$qPtr#9105][368LL])) then
                particles_pushElement($qBasePtr12#9038, $idx#9104, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q12#8846[$qPtr#9105][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9104 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9106 : int3d = int3d({1, 1, -1})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9106)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9107 : int32 = 0
            for $qPtr#9108 : int1d(int8[376], $q13#8858) in $q13#8858 do
              if (not bool($q13#8858[$qPtr#9108][368LL])) then
                particles_pushElement($qBasePtr13#9040, $idx#9107, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q13#8858[$qPtr#9108][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9107 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9109 : int3d = int3d({1, -1, 0})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9109)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9110 : int32 = 0
            for $qPtr#9111 : int1d(int8[376], $q14#8870) in $q14#8870 do
              if (not bool($q14#8870[$qPtr#9111][368LL])) then
                particles_pushElement($qBasePtr14#9042, $idx#9110, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q14#8870[$qPtr#9111][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9110 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9112 : int3d = int3d({1, -1, 1})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9112)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9113 : int32 = 0
            for $qPtr#9114 : int1d(int8[376], $q15#8882) in $q15#8882 do
              if (not bool($q15#8882[$qPtr#9114][368LL])) then
                particles_pushElement($qBasePtr15#9044, $idx#9113, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q15#8882[$qPtr#9114][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9113 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9115 : int3d = int3d({1, -1, -1})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9115)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9116 : int32 = 0
            for $qPtr#9117 : int1d(int8[376], $q16#8894) in $q16#8894 do
              if (not bool($q16#8894[$qPtr#9117][368LL])) then
                particles_pushElement($qBasePtr16#9046, $idx#9116, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q16#8894[$qPtr#9117][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9116 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9118 : int3d = int3d({-1, 0, 0})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9118)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9119 : int32 = 0
            for $qPtr#9120 : int1d(int8[376], $q17#8906) in $q17#8906 do
              if (not bool($q17#8906[$qPtr#9120][368LL])) then
                particles_pushElement($qBasePtr17#9048, $idx#9119, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q17#8906[$qPtr#9120][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9119 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9121 : int3d = int3d({-1, 0, 1})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9121)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9122 : int32 = 0
            for $qPtr#9123 : int1d(int8[376], $q18#8918) in $q18#8918 do
              if (not bool($q18#8918[$qPtr#9123][368LL])) then
                particles_pushElement($qBasePtr18#9050, $idx#9122, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q18#8918[$qPtr#9123][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9122 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9124 : int3d = int3d({-1, 0, -1})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9124)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9125 : int32 = 0
            for $qPtr#9126 : int1d(int8[376], $q19#8930) in $q19#8930 do
              if (not bool($q19#8930[$qPtr#9126][368LL])) then
                particles_pushElement($qBasePtr19#9052, $idx#9125, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q19#8930[$qPtr#9126][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9125 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9127 : int3d = int3d({-1, 1, 0})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9127)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9128 : int32 = 0
            for $qPtr#9129 : int1d(int8[376], $q20#8942) in $q20#8942 do
              if (not bool($q20#8942[$qPtr#9129][368LL])) then
                particles_pushElement($qBasePtr20#9054, $idx#9128, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q20#8942[$qPtr#9129][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9128 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9130 : int3d = int3d({-1, 1, 1})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9130)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9131 : int32 = 0
            for $qPtr#9132 : int1d(int8[376], $q21#8954) in $q21#8954 do
              if (not bool($q21#8954[$qPtr#9132][368LL])) then
                particles_pushElement($qBasePtr21#9056, $idx#9131, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q21#8954[$qPtr#9132][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9131 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9133 : int3d = int3d({-1, 1, -1})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9133)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9134 : int32 = 0
            for $qPtr#9135 : int1d(int8[376], $q22#8966) in $q22#8966 do
              if (not bool($q22#8966[$qPtr#9135][368LL])) then
                particles_pushElement($qBasePtr22#9058, $idx#9134, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q22#8966[$qPtr#9135][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9134 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9136 : int3d = int3d({-1, -1, 0})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9136)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9137 : int32 = 0
            for $qPtr#9138 : int1d(int8[376], $q23#8978) in $q23#8978 do
              if (not bool($q23#8978[$qPtr#9138][368LL])) then
                particles_pushElement($qBasePtr23#9060, $idx#9137, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q23#8978[$qPtr#9138][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9137 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9139 : int3d = int3d({-1, -1, 1})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9139)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9140 : int32 = 0
            for $qPtr#9141 : int1d(int8[376], $q24#8990) in $q24#8990 do
              if (not bool($q24#8990[$qPtr#9141][368LL])) then
                particles_pushElement($qBasePtr24#9062, $idx#9140, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q24#8990[$qPtr#9141][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9140 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        do
          var $colorOff#9142 : int3d = int3d({-1, -1, -1})
          if ((@$rPtr#9065).__valid and ($elemColor#9066==((($partColor#8698+$colorOff#9142)+{2, 3, 4})%{2, 3, 4}))) then
            var $idx#9143 : int32 = 0
            for $qPtr#9144 : int1d(int8[376], $q25#9002) in $q25#9002 do
              if (not bool($q25#9002[$qPtr#9144][368LL])) then
                particles_pushElement($qBasePtr25#9064, $idx#9143, $r#8696[$rPtr#9065])
                (@$rPtr#9065).__valid = false
                std.assert(bool($q25#9002[$qPtr#9144][368LL]), "Element did not get copied properly")
                break
              else
              end
              $idx#9143 += 1
            end
            std.assert((not (@$rPtr#9065).__valid), "Transfer queue ran out of space")
          else
          end
        end
        std.assert((not (@$rPtr#9065).__valid), "Element moved past predicted stencil")
      else
      end
    else
    end
  end
end
terra particles_pullElement(src : &int8) : particles_columns
    var dst : particles_columns
    memcpy([&opaque](&dst), [&opaque](src), [uint64](376))
    return dst
end
task particles_pullAll($color#10511 : int3d, $r#10512 : region#2667(ispace#2667(int1d), particles_columns), $q0#10457 : region#2638(ispace#2638(int1d), int8[376]), $q1#10459 : region#2639(ispace#2639(int1d), int8[376]), $q2#10461 : region#2640(ispace#2640(int1d), int8[376]), $q3#10463 : region#2641(ispace#2641(int1d), int8[376]), $q4#10465 : region#2642(ispace#2642(int1d), int8[376]), $q5#10467 : region#2643(ispace#2643(int1d), int8[376]), $q6#10469 : region#2644(ispace#2644(int1d), int8[376]), $q7#10471 : region#2645(ispace#2645(int1d), int8[376]), $q8#10473 : region#2646(ispace#2646(int1d), int8[376]), $q9#10475 : region#2647(ispace#2647(int1d), int8[376]), $q10#10477 : region#2648(ispace#2648(int1d), int8[376]), $q11#10479 : region#2649(ispace#2649(int1d), int8[376]), $q12#10481 : region#2650(ispace#2650(int1d), int8[376]), $q13#10483 : region#2651(ispace#2651(int1d), int8[376]), $q14#10485 : region#2652(ispace#2652(int1d), int8[376]), $q15#10487 : region#2653(ispace#2653(int1d), int8[376]), $q16#10489 : region#2654(ispace#2654(int1d), int8[376]), $q17#10491 : region#2655(ispace#2655(int1d), int8[376]), $q18#10493 : region#2656(ispace#2656(int1d), int8[376]), $q19#10495 : region#2657(ispace#2657(int1d), int8[376]), $q20#10497 : region#2658(ispace#2658(int1d), int8[376]), $q21#10499 : region#2659(ispace#2659(int1d), int8[376]), $q22#10501 : region#2660(ispace#2660(int1d), int8[376]), $q23#10503 : region#2661(ispace#2661(int1d), int8[376]), $q24#10505 : region#2662(ispace#2662(int1d), int8[376]), $q25#10507 : region#2663(ispace#2663(int1d), int8[376]))
-- leaf (false), inner (false), idempotent (false)
where
  reads($r#10512), writes($r#10512), reads($q0#10457), reads($q1#10459), reads($q2#10461), reads($q3#10463), reads($q4#10465), reads($q5#10467), reads($q6#10469), reads($q7#10471), reads($q8#10473), reads($q9#10475), reads($q10#10477), reads($q11#10479), reads($q12#10481), reads($q13#10483), reads($q14#10485), reads($q15#10487), reads($q16#10489), reads($q17#10491), reads($q18#10493), reads($q19#10495), reads($q20#10497), reads($q21#10499), reads($q22#10501), reads($q23#10503), reads($q24#10505), reads($q25#10507)
do
  for $qPtr#10670 : int1d(int8[376], $q0#10457) in $q0#10457 do
    if bool($q0#10457[$qPtr#10670][368LL]) then
      var $copied#10671 : bool = false
      for $rPtr#10672 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10672).__valid) then
          $r#10512[$rPtr#10672] = particles_pullElement(&int8($q0#10457[$qPtr#10670]))
          $copied#10671 = true
          std.assert($r#10512[$rPtr#10672].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10671, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10673 : int1d(int8[376], $q1#10459) in $q1#10459 do
    if bool($q1#10459[$qPtr#10673][368LL]) then
      var $copied#10674 : bool = false
      for $rPtr#10675 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10675).__valid) then
          $r#10512[$rPtr#10675] = particles_pullElement(&int8($q1#10459[$qPtr#10673]))
          $copied#10674 = true
          std.assert($r#10512[$rPtr#10675].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10674, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10676 : int1d(int8[376], $q2#10461) in $q2#10461 do
    if bool($q2#10461[$qPtr#10676][368LL]) then
      var $copied#10677 : bool = false
      for $rPtr#10678 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10678).__valid) then
          $r#10512[$rPtr#10678] = particles_pullElement(&int8($q2#10461[$qPtr#10676]))
          $copied#10677 = true
          std.assert($r#10512[$rPtr#10678].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10677, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10679 : int1d(int8[376], $q3#10463) in $q3#10463 do
    if bool($q3#10463[$qPtr#10679][368LL]) then
      var $copied#10680 : bool = false
      for $rPtr#10681 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10681).__valid) then
          $r#10512[$rPtr#10681] = particles_pullElement(&int8($q3#10463[$qPtr#10679]))
          $copied#10680 = true
          std.assert($r#10512[$rPtr#10681].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10680, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10682 : int1d(int8[376], $q4#10465) in $q4#10465 do
    if bool($q4#10465[$qPtr#10682][368LL]) then
      var $copied#10683 : bool = false
      for $rPtr#10684 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10684).__valid) then
          $r#10512[$rPtr#10684] = particles_pullElement(&int8($q4#10465[$qPtr#10682]))
          $copied#10683 = true
          std.assert($r#10512[$rPtr#10684].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10683, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10685 : int1d(int8[376], $q5#10467) in $q5#10467 do
    if bool($q5#10467[$qPtr#10685][368LL]) then
      var $copied#10686 : bool = false
      for $rPtr#10687 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10687).__valid) then
          $r#10512[$rPtr#10687] = particles_pullElement(&int8($q5#10467[$qPtr#10685]))
          $copied#10686 = true
          std.assert($r#10512[$rPtr#10687].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10686, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10688 : int1d(int8[376], $q6#10469) in $q6#10469 do
    if bool($q6#10469[$qPtr#10688][368LL]) then
      var $copied#10689 : bool = false
      for $rPtr#10690 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10690).__valid) then
          $r#10512[$rPtr#10690] = particles_pullElement(&int8($q6#10469[$qPtr#10688]))
          $copied#10689 = true
          std.assert($r#10512[$rPtr#10690].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10689, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10691 : int1d(int8[376], $q7#10471) in $q7#10471 do
    if bool($q7#10471[$qPtr#10691][368LL]) then
      var $copied#10692 : bool = false
      for $rPtr#10693 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10693).__valid) then
          $r#10512[$rPtr#10693] = particles_pullElement(&int8($q7#10471[$qPtr#10691]))
          $copied#10692 = true
          std.assert($r#10512[$rPtr#10693].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10692, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10694 : int1d(int8[376], $q8#10473) in $q8#10473 do
    if bool($q8#10473[$qPtr#10694][368LL]) then
      var $copied#10695 : bool = false
      for $rPtr#10696 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10696).__valid) then
          $r#10512[$rPtr#10696] = particles_pullElement(&int8($q8#10473[$qPtr#10694]))
          $copied#10695 = true
          std.assert($r#10512[$rPtr#10696].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10695, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10697 : int1d(int8[376], $q9#10475) in $q9#10475 do
    if bool($q9#10475[$qPtr#10697][368LL]) then
      var $copied#10698 : bool = false
      for $rPtr#10699 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10699).__valid) then
          $r#10512[$rPtr#10699] = particles_pullElement(&int8($q9#10475[$qPtr#10697]))
          $copied#10698 = true
          std.assert($r#10512[$rPtr#10699].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10698, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10700 : int1d(int8[376], $q10#10477) in $q10#10477 do
    if bool($q10#10477[$qPtr#10700][368LL]) then
      var $copied#10701 : bool = false
      for $rPtr#10702 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10702).__valid) then
          $r#10512[$rPtr#10702] = particles_pullElement(&int8($q10#10477[$qPtr#10700]))
          $copied#10701 = true
          std.assert($r#10512[$rPtr#10702].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10701, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10703 : int1d(int8[376], $q11#10479) in $q11#10479 do
    if bool($q11#10479[$qPtr#10703][368LL]) then
      var $copied#10704 : bool = false
      for $rPtr#10705 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10705).__valid) then
          $r#10512[$rPtr#10705] = particles_pullElement(&int8($q11#10479[$qPtr#10703]))
          $copied#10704 = true
          std.assert($r#10512[$rPtr#10705].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10704, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10706 : int1d(int8[376], $q12#10481) in $q12#10481 do
    if bool($q12#10481[$qPtr#10706][368LL]) then
      var $copied#10707 : bool = false
      for $rPtr#10708 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10708).__valid) then
          $r#10512[$rPtr#10708] = particles_pullElement(&int8($q12#10481[$qPtr#10706]))
          $copied#10707 = true
          std.assert($r#10512[$rPtr#10708].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10707, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10709 : int1d(int8[376], $q13#10483) in $q13#10483 do
    if bool($q13#10483[$qPtr#10709][368LL]) then
      var $copied#10710 : bool = false
      for $rPtr#10711 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10711).__valid) then
          $r#10512[$rPtr#10711] = particles_pullElement(&int8($q13#10483[$qPtr#10709]))
          $copied#10710 = true
          std.assert($r#10512[$rPtr#10711].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10710, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10712 : int1d(int8[376], $q14#10485) in $q14#10485 do
    if bool($q14#10485[$qPtr#10712][368LL]) then
      var $copied#10713 : bool = false
      for $rPtr#10714 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10714).__valid) then
          $r#10512[$rPtr#10714] = particles_pullElement(&int8($q14#10485[$qPtr#10712]))
          $copied#10713 = true
          std.assert($r#10512[$rPtr#10714].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10713, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10715 : int1d(int8[376], $q15#10487) in $q15#10487 do
    if bool($q15#10487[$qPtr#10715][368LL]) then
      var $copied#10716 : bool = false
      for $rPtr#10717 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10717).__valid) then
          $r#10512[$rPtr#10717] = particles_pullElement(&int8($q15#10487[$qPtr#10715]))
          $copied#10716 = true
          std.assert($r#10512[$rPtr#10717].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10716, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10718 : int1d(int8[376], $q16#10489) in $q16#10489 do
    if bool($q16#10489[$qPtr#10718][368LL]) then
      var $copied#10719 : bool = false
      for $rPtr#10720 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10720).__valid) then
          $r#10512[$rPtr#10720] = particles_pullElement(&int8($q16#10489[$qPtr#10718]))
          $copied#10719 = true
          std.assert($r#10512[$rPtr#10720].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10719, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10721 : int1d(int8[376], $q17#10491) in $q17#10491 do
    if bool($q17#10491[$qPtr#10721][368LL]) then
      var $copied#10722 : bool = false
      for $rPtr#10723 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10723).__valid) then
          $r#10512[$rPtr#10723] = particles_pullElement(&int8($q17#10491[$qPtr#10721]))
          $copied#10722 = true
          std.assert($r#10512[$rPtr#10723].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10722, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10724 : int1d(int8[376], $q18#10493) in $q18#10493 do
    if bool($q18#10493[$qPtr#10724][368LL]) then
      var $copied#10725 : bool = false
      for $rPtr#10726 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10726).__valid) then
          $r#10512[$rPtr#10726] = particles_pullElement(&int8($q18#10493[$qPtr#10724]))
          $copied#10725 = true
          std.assert($r#10512[$rPtr#10726].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10725, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10727 : int1d(int8[376], $q19#10495) in $q19#10495 do
    if bool($q19#10495[$qPtr#10727][368LL]) then
      var $copied#10728 : bool = false
      for $rPtr#10729 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10729).__valid) then
          $r#10512[$rPtr#10729] = particles_pullElement(&int8($q19#10495[$qPtr#10727]))
          $copied#10728 = true
          std.assert($r#10512[$rPtr#10729].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10728, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10730 : int1d(int8[376], $q20#10497) in $q20#10497 do
    if bool($q20#10497[$qPtr#10730][368LL]) then
      var $copied#10731 : bool = false
      for $rPtr#10732 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10732).__valid) then
          $r#10512[$rPtr#10732] = particles_pullElement(&int8($q20#10497[$qPtr#10730]))
          $copied#10731 = true
          std.assert($r#10512[$rPtr#10732].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10731, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10733 : int1d(int8[376], $q21#10499) in $q21#10499 do
    if bool($q21#10499[$qPtr#10733][368LL]) then
      var $copied#10734 : bool = false
      for $rPtr#10735 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10735).__valid) then
          $r#10512[$rPtr#10735] = particles_pullElement(&int8($q21#10499[$qPtr#10733]))
          $copied#10734 = true
          std.assert($r#10512[$rPtr#10735].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10734, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10736 : int1d(int8[376], $q22#10501) in $q22#10501 do
    if bool($q22#10501[$qPtr#10736][368LL]) then
      var $copied#10737 : bool = false
      for $rPtr#10738 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10738).__valid) then
          $r#10512[$rPtr#10738] = particles_pullElement(&int8($q22#10501[$qPtr#10736]))
          $copied#10737 = true
          std.assert($r#10512[$rPtr#10738].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10737, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10739 : int1d(int8[376], $q23#10503) in $q23#10503 do
    if bool($q23#10503[$qPtr#10739][368LL]) then
      var $copied#10740 : bool = false
      for $rPtr#10741 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10741).__valid) then
          $r#10512[$rPtr#10741] = particles_pullElement(&int8($q23#10503[$qPtr#10739]))
          $copied#10740 = true
          std.assert($r#10512[$rPtr#10741].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10740, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10742 : int1d(int8[376], $q24#10505) in $q24#10505 do
    if bool($q24#10505[$qPtr#10742][368LL]) then
      var $copied#10743 : bool = false
      for $rPtr#10744 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10744).__valid) then
          $r#10512[$rPtr#10744] = particles_pullElement(&int8($q24#10505[$qPtr#10742]))
          $copied#10743 = true
          std.assert($r#10512[$rPtr#10744].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10743, "Ran out of space on sub-partition")
    else
    end
  end
  for $qPtr#10745 : int1d(int8[376], $q25#10507) in $q25#10507 do
    if bool($q25#10507[$qPtr#10745][368LL]) then
      var $copied#10746 : bool = false
      for $rPtr#10747 : int1d(particles_columns, $r#10512) in $r#10512 do
        if (not (@$rPtr#10747).__valid) then
          $r#10512[$rPtr#10747] = particles_pullElement(&int8($q25#10507[$qPtr#10745]))
          $copied#10746 = true
          std.assert($r#10512[$rPtr#10747].__valid, "Pulled particle was not copied correctly")
          break
        else
        end
      end
      std.assert($copied#10746, "Ran out of space on sub-partition")
    else
    end
  end
end
task TrilinearInterpolateVelocity($xyz#11487 : double[3], $c000#11488 : double[3], $c100#11489 : double[3], $c010#11490 : double[3], $c110#11491 : double[3], $c001#11492 : double[3], $c101#11493 : double[3], $c011#11494 : double[3], $c111#11495 : double[3]) : double[3]
-- leaf (false), inner (false), idempotent (false)
  var $dX#11532 : double = fmod(((($xyz#11487[int32(0)]-double(0.6987052))/double(0.0015841600529101))+double(0.5)), double(1))
  var $dY#11533 : double = fmod(((($xyz#11487[int32(1)]-double(0.99700342743682))/double(0.0010188725631769))+double(0.5)), double(1))
  var $dZ#11534 : double = fmod(((($xyz#11487[int32(2)]-double(0.70456500505464))/double(0.0011352949453552))+double(0.5)), double(1))
  var $oneMinusdX#11535 : double = (double(1)-$dX#11532)
  var $oneMinusdY#11536 : double = (double(1)-$dY#11533)
  var $oneMinusdZ#11537 : double = (double(1)-$dZ#11534)
  var $weight00#11538 : double[3] = vv_add_double_3(vs_mul_double_3($c000#11488, $oneMinusdX#11535), vs_mul_double_3($c100#11489, $dX#11532))
  var $weight10#11539 : double[3] = vv_add_double_3(vs_mul_double_3($c010#11490, $oneMinusdX#11535), vs_mul_double_3($c110#11491, $dX#11532))
  var $weight01#11540 : double[3] = vv_add_double_3(vs_mul_double_3($c001#11492, $oneMinusdX#11535), vs_mul_double_3($c101#11493, $dX#11532))
  var $weight11#11541 : double[3] = vv_add_double_3(vs_mul_double_3($c011#11494, $oneMinusdX#11535), vs_mul_double_3($c111#11495, $dX#11532))
  var $weight0#11542 : double[3] = vv_add_double_3(vs_mul_double_3($weight00#11538, $oneMinusdY#11536), vs_mul_double_3($weight10#11539, $dY#11533))
  var $weight1#11543 : double[3] = vv_add_double_3(vs_mul_double_3($weight01#11540, $oneMinusdY#11536), vs_mul_double_3($weight11#11541, $dY#11533))
  return vv_add_double_3(vs_mul_double_3($weight0#11542, $oneMinusdZ#11537), vs_mul_double_3($weight1#11543, $dZ#11534))
end
task InterpolateTriVelocity($c#11390 : int3d, $xyz#11391 : double[3], $Fluid#11393 : region#2855(ispace#2855(int3d), Fluid_columns)) : double[3]
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#11393.centerCoordinates), reads($Fluid#11393.velocity)
do
  var $velocity000#11645 : double[3] = double[3](array(double(0), double(0), double(0)))
  var $velocity100#11646 : double[3] = double[3](array(double(0), double(0), double(0)))
  var $velocity010#11647 : double[3] = double[3](array(double(0), double(0), double(0)))
  var $velocity110#11648 : double[3] = double[3](array(double(0), double(0), double(0)))
  var $velocity001#11649 : double[3] = double[3](array(double(0), double(0), double(0)))
  var $velocity101#11650 : double[3] = double[3](array(double(0), double(0), double(0)))
  var $velocity011#11651 : double[3] = double[3](array(double(0), double(0), double(0)))
  var $velocity111#11652 : double[3] = double[3](array(double(0), double(0), double(0)))
  var $velocity0#11653 : double[3] = $Fluid#11393[$c#11390].velocity
  if ($xyz#11391[int32(0)]>$Fluid#11393[$c#11390].centerCoordinates[int32(0)]) then
    var $velocityb#11654 : double[3] = $Fluid#11393[(($c#11390+{1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
    if ($xyz#11391[int32(1)]>$Fluid#11393[$c#11390].centerCoordinates[int32(1)]) then
      var $velocityaa#11655 : double[3] = $Fluid#11393[(($c#11390+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
      var $velocityab#11656 : double[3] = $Fluid#11393[(($c#11390+{1, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
      if ($xyz#11391[int32(2)]>$Fluid#11393[$c#11390].centerCoordinates[int32(2)]) then
        $velocity000#11645 = $velocity0#11653
        $velocity100#11646 = $velocityb#11654
        $velocity010#11647 = $velocityaa#11655
        $velocity110#11648 = $velocityab#11656
        $velocity001#11649 = $Fluid#11393[(($c#11390+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity101#11650 = $Fluid#11393[(($c#11390+{1, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity011#11651 = $Fluid#11393[(($c#11390+{0, 1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity111#11652 = $Fluid#11393[(($c#11390+{1, 1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
      else
        $velocity000#11645 = $Fluid#11393[(($c#11390+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity100#11646 = $Fluid#11393[(($c#11390+{1, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity010#11647 = $Fluid#11393[(($c#11390+{0, 1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity110#11648 = $Fluid#11393[(($c#11390+{1, 1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity001#11649 = $velocity0#11653
        $velocity101#11650 = $velocityb#11654
        $velocity011#11651 = $velocityaa#11655
        $velocity111#11652 = $velocityab#11656
      end
    else
      var $velocityaa#11657 : double[3] = $Fluid#11393[(($c#11390+{0, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
      var $velocityab#11658 : double[3] = $Fluid#11393[(($c#11390+{1, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
      if ($xyz#11391[int32(2)]>$Fluid#11393[$c#11390].centerCoordinates[int32(2)]) then
        $velocity000#11645 = $velocityaa#11657
        $velocity100#11646 = $velocityab#11658
        $velocity010#11647 = $velocity0#11653
        $velocity110#11648 = $velocityb#11654
        $velocity001#11649 = $Fluid#11393[(($c#11390+{0, -1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity101#11650 = $Fluid#11393[(($c#11390+{1, -1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity011#11651 = $Fluid#11393[(($c#11390+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity111#11652 = $Fluid#11393[(($c#11390+{1, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
      else
        $velocity000#11645 = $Fluid#11393[(($c#11390+{0, -1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity100#11646 = $Fluid#11393[(($c#11390+{1, -1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity010#11647 = $Fluid#11393[(($c#11390+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity110#11648 = $Fluid#11393[(($c#11390+{1, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity001#11649 = $velocityaa#11657
        $velocity101#11650 = $velocityab#11658
        $velocity011#11651 = $velocity0#11653
        $velocity111#11652 = $velocityb#11654
      end
    end
  else
    var $velocitya#11659 : double[3] = $Fluid#11393[(($c#11390+{-1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
    if ($xyz#11391[int32(1)]>$Fluid#11393[$c#11390].centerCoordinates[int32(1)]) then
      var $velocityaa#11660 : double[3] = $Fluid#11393[(($c#11390+{-1, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
      var $velocityab#11661 : double[3] = $Fluid#11393[(($c#11390+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
      if ($xyz#11391[int32(2)]>$Fluid#11393[$c#11390].centerCoordinates[int32(2)]) then
        $velocity000#11645 = $velocitya#11659
        $velocity100#11646 = $velocity0#11653
        $velocity010#11647 = $velocityaa#11660
        $velocity110#11648 = $velocityab#11661
        $velocity001#11649 = $Fluid#11393[(($c#11390+{-1, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity101#11650 = $Fluid#11393[(($c#11390+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity011#11651 = $Fluid#11393[(($c#11390+{-1, 1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity111#11652 = $Fluid#11393[(($c#11390+{0, 1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
      else
        $velocity000#11645 = $Fluid#11393[(($c#11390+{-1, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity100#11646 = $Fluid#11393[(($c#11390+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity010#11647 = $Fluid#11393[(($c#11390+{-1, 1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity110#11648 = $Fluid#11393[(($c#11390+{0, 1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity001#11649 = $velocitya#11659
        $velocity101#11650 = $velocity0#11653
        $velocity011#11651 = $velocityaa#11660
        $velocity111#11652 = $velocityab#11661
      end
    else
      var $velocityaa#11662 : double[3] = $Fluid#11393[(($c#11390+{-1, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
      var $velocityab#11663 : double[3] = $Fluid#11393[(($c#11390+{0, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
      if ($xyz#11391[int32(2)]>$Fluid#11393[$c#11390].centerCoordinates[int32(2)]) then
        $velocity000#11645 = $velocityaa#11662
        $velocity100#11646 = $velocityab#11663
        $velocity010#11647 = $velocitya#11659
        $velocity110#11648 = $velocity0#11653
        $velocity001#11649 = $Fluid#11393[(($c#11390+{-1, -1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity101#11650 = $Fluid#11393[(($c#11390+{0, -1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity011#11651 = $Fluid#11393[(($c#11390+{-1, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity111#11652 = $Fluid#11393[(($c#11390+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
      else
        $velocity000#11645 = $Fluid#11393[(($c#11390+{-1, -1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity100#11646 = $Fluid#11393[(($c#11390+{0, -1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity010#11647 = $Fluid#11393[(($c#11390+{-1, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity110#11648 = $Fluid#11393[(($c#11390+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
        $velocity001#11649 = $velocityaa#11662
        $velocity101#11650 = $velocityab#11663
        $velocity011#11651 = $velocitya#11659
        $velocity111#11652 = $velocity0#11653
      end
    end
  end
  var $11740 : double[3]
  do
    var $11731 : double[3] = $xyz#11391
    var $11732 : double[3] = $velocity000#11645
    var $11733 : double[3] = $velocity100#11646
    var $11734 : double[3] = $velocity010#11647
    var $11735 : double[3] = $velocity110#11648
    var $11736 : double[3] = $velocity001#11649
    var $11737 : double[3] = $velocity101#11650
    var $11738 : double[3] = $velocity011#11651
    var $11739 : double[3] = $velocity111#11652
    var $dX#11532 : double = fmod(((($11731[int32(0)]-double(0.6987052))/double(0.0015841600529101))+double(0.5)), double(1))
    var $dY#11533 : double = fmod(((($11731[int32(1)]-double(0.99700342743682))/double(0.0010188725631769))+double(0.5)), double(1))
    var $dZ#11534 : double = fmod(((($11731[int32(2)]-double(0.70456500505464))/double(0.0011352949453552))+double(0.5)), double(1))
    var $oneMinusdX#11535 : double = (double(1)-$dX#11532)
    var $oneMinusdY#11536 : double = (double(1)-$dY#11533)
    var $oneMinusdZ#11537 : double = (double(1)-$dZ#11534)
    var $weight00#11538 : double[3] = vv_add_double_3(vs_mul_double_3($11732, $oneMinusdX#11535), vs_mul_double_3($11733, $dX#11532))
    var $weight10#11539 : double[3] = vv_add_double_3(vs_mul_double_3($11734, $oneMinusdX#11535), vs_mul_double_3($11735, $dX#11532))
    var $weight01#11540 : double[3] = vv_add_double_3(vs_mul_double_3($11736, $oneMinusdX#11535), vs_mul_double_3($11737, $dX#11532))
    var $weight11#11541 : double[3] = vv_add_double_3(vs_mul_double_3($11738, $oneMinusdX#11535), vs_mul_double_3($11739, $dX#11532))
    var $weight0#11542 : double[3] = vv_add_double_3(vs_mul_double_3($weight00#11538, $oneMinusdY#11536), vs_mul_double_3($weight10#11539, $dY#11533))
    var $weight1#11543 : double[3] = vv_add_double_3(vs_mul_double_3($weight01#11540, $oneMinusdY#11536), vs_mul_double_3($weight11#11541, $dY#11533))
    $11740 = vv_add_double_3(vs_mul_double_3($weight0#11542, $oneMinusdZ#11537), vs_mul_double_3($weight1#11543, $dZ#11534))
  end
  var $dummy#11741 : int32 = 0
  return $11740
end
task TrilinearInterpolateTemp($xyz#11894 : double[3], $c000#11895 : double, $c100#11896 : double, $c010#11897 : double, $c110#11898 : double, $c001#11899 : double, $c101#11900 : double, $c011#11901 : double, $c111#11902 : double) : double
-- leaf (false), inner (false), idempotent (false)
  var $dX#11939 : double = fmod(((($xyz#11894[int32(0)]-double(0.6987052))/double(0.0015841600529101))+double(0.5)), double(1))
  var $dY#11940 : double = fmod(((($xyz#11894[int32(1)]-double(0.99700342743682))/double(0.0010188725631769))+double(0.5)), double(1))
  var $dZ#11941 : double = fmod(((($xyz#11894[int32(2)]-double(0.70456500505464))/double(0.0011352949453552))+double(0.5)), double(1))
  var $oneMinusdX#11942 : double = (double(1)-$dX#11939)
  var $oneMinusdY#11943 : double = (double(1)-$dY#11940)
  var $oneMinusdZ#11944 : double = (double(1)-$dZ#11941)
  var $weight00#11945 : double = (($c000#11895*$oneMinusdX#11942)+($c100#11896*$dX#11939))
  var $weight10#11946 : double = (($c010#11897*$oneMinusdX#11942)+($c110#11898*$dX#11939))
  var $weight01#11947 : double = (($c001#11899*$oneMinusdX#11942)+($c101#11900*$dX#11939))
  var $weight11#11948 : double = (($c011#11901*$oneMinusdX#11942)+($c111#11902*$dX#11939))
  var $weight0#11949 : double = (($weight00#11945*$oneMinusdY#11943)+($weight10#11946*$dY#11940))
  var $weight1#11950 : double = (($weight01#11947*$oneMinusdY#11943)+($weight11#11948*$dY#11940))
  return (($weight0#11949*$oneMinusdZ#11944)+($weight1#11950*$dZ#11941))
end
task InterpolateTriTemp($c#11797 : int3d, $xyz#11798 : double[3], $Fluid#11800 : region#2925(ispace#2925(int3d), Fluid_columns)) : double
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#11800.centerCoordinates), reads($Fluid#11800.temperature)
do
  var $temp000#11996 : double = double(double(0))
  var $temp100#11997 : double = double(double(0))
  var $temp010#11998 : double = double(double(0))
  var $temp110#11999 : double = double(double(0))
  var $temp001#12000 : double = double(double(0))
  var $temp101#12001 : double = double(double(0))
  var $temp011#12002 : double = double(double(0))
  var $temp111#12003 : double = double(double(0))
  var $temp0#12004 : double = $Fluid#11800[$c#11797].temperature
  if ($xyz#11798[int32(0)]>$Fluid#11800[$c#11797].centerCoordinates[int32(0)]) then
    var $tempb#12005 : double = $Fluid#11800[(($c#11797+{1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
    if ($xyz#11798[int32(1)]>$Fluid#11800[$c#11797].centerCoordinates[int32(1)]) then
      var $tempaa#12006 : double = $Fluid#11800[(($c#11797+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
      var $tempab#12007 : double = $Fluid#11800[(($c#11797+{1, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
      if ($xyz#11798[int32(2)]>$Fluid#11800[$c#11797].centerCoordinates[int32(2)]) then
        $temp000#11996 = $temp0#12004
        $temp100#11997 = $tempb#12005
        $temp010#11998 = $tempaa#12006
        $temp110#11999 = $tempab#12007
        $temp001#12000 = $Fluid#11800[(($c#11797+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp101#12001 = $Fluid#11800[(($c#11797+{1, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp011#12002 = $Fluid#11800[(($c#11797+{0, 1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp111#12003 = $Fluid#11800[(($c#11797+{1, 1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
      else
        $temp000#11996 = $Fluid#11800[(($c#11797+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp100#11997 = $Fluid#11800[(($c#11797+{1, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp010#11998 = $Fluid#11800[(($c#11797+{0, 1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp110#11999 = $Fluid#11800[(($c#11797+{1, 1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp001#12000 = $temp0#12004
        $temp101#12001 = $tempb#12005
        $temp011#12002 = $tempaa#12006
        $temp111#12003 = $tempab#12007
      end
    else
      var $tempaa#12008 : double = $Fluid#11800[(($c#11797+{0, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
      var $tempab#12009 : double = $Fluid#11800[(($c#11797+{1, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
      if ($xyz#11798[int32(2)]>$Fluid#11800[$c#11797].centerCoordinates[int32(2)]) then
        $temp000#11996 = $tempaa#12008
        $temp100#11997 = $tempab#12009
        $temp010#11998 = $temp0#12004
        $temp110#11999 = $tempb#12005
        $temp001#12000 = $Fluid#11800[(($c#11797+{0, -1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp101#12001 = $Fluid#11800[(($c#11797+{1, -1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp011#12002 = $Fluid#11800[(($c#11797+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp111#12003 = $Fluid#11800[(($c#11797+{1, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
      else
        $temp000#11996 = $Fluid#11800[(($c#11797+{0, -1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp100#11997 = $Fluid#11800[(($c#11797+{1, -1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp010#11998 = $Fluid#11800[(($c#11797+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp110#11999 = $Fluid#11800[(($c#11797+{1, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp001#12000 = $tempaa#12008
        $temp101#12001 = $tempab#12009
        $temp011#12002 = $temp0#12004
        $temp111#12003 = $tempb#12005
      end
    end
  else
    var $tempa#12010 : double = $Fluid#11800[(($c#11797+{-1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
    if ($xyz#11798[int32(1)]>$Fluid#11800[$c#11797].centerCoordinates[int32(1)]) then
      var $tempaa#12011 : double = $Fluid#11800[(($c#11797+{-1, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
      var $tempab#12012 : double = $Fluid#11800[(($c#11797+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
      if ($xyz#11798[int32(2)]>$Fluid#11800[$c#11797].centerCoordinates[int32(2)]) then
        $temp000#11996 = $tempa#12010
        $temp100#11997 = $temp0#12004
        $temp010#11998 = $tempaa#12011
        $temp110#11999 = $tempab#12012
        $temp001#12000 = $Fluid#11800[(($c#11797+{-1, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp101#12001 = $Fluid#11800[(($c#11797+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp011#12002 = $Fluid#11800[(($c#11797+{-1, 1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp111#12003 = $Fluid#11800[(($c#11797+{0, 1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
      else
        $temp000#11996 = $Fluid#11800[(($c#11797+{-1, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp100#11997 = $Fluid#11800[(($c#11797+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp010#11998 = $Fluid#11800[(($c#11797+{-1, 1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp110#11999 = $Fluid#11800[(($c#11797+{0, 1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp001#12000 = $tempa#12010
        $temp101#12001 = $temp0#12004
        $temp011#12002 = $tempaa#12011
        $temp111#12003 = $tempab#12012
      end
    else
      var $tempaa#12013 : double = $Fluid#11800[(($c#11797+{-1, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
      var $tempab#12014 : double = $Fluid#11800[(($c#11797+{0, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
      if ($xyz#11798[int32(2)]>$Fluid#11800[$c#11797].centerCoordinates[int32(2)]) then
        $temp000#11996 = $tempaa#12013
        $temp100#11997 = $tempab#12014
        $temp010#11998 = $tempa#12010
        $temp110#11999 = $temp0#12004
        $temp001#12000 = $Fluid#11800[(($c#11797+{-1, -1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp101#12001 = $Fluid#11800[(($c#11797+{0, -1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp011#12002 = $Fluid#11800[(($c#11797+{-1, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp111#12003 = $Fluid#11800[(($c#11797+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
      else
        $temp000#11996 = $Fluid#11800[(($c#11797+{-1, -1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp100#11997 = $Fluid#11800[(($c#11797+{0, -1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp010#11998 = $Fluid#11800[(($c#11797+{-1, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp110#11999 = $Fluid#11800[(($c#11797+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
        $temp001#12000 = $tempaa#12013
        $temp101#12001 = $tempab#12014
        $temp011#12002 = $tempa#12010
        $temp111#12003 = $temp0#12004
      end
    end
  end
  var $12043 : double
  do
    var $12034 : double[3] = $xyz#11798
    var $12035 : double = $temp000#11996
    var $12036 : double = $temp100#11997
    var $12037 : double = $temp010#11998
    var $12038 : double = $temp110#11999
    var $12039 : double = $temp001#12000
    var $12040 : double = $temp101#12001
    var $12041 : double = $temp011#12002
    var $12042 : double = $temp111#12003
    var $dX#11939 : double = fmod(((($12034[int32(0)]-double(0.6987052))/double(0.0015841600529101))+double(0.5)), double(1))
    var $dY#11940 : double = fmod(((($12034[int32(1)]-double(0.99700342743682))/double(0.0010188725631769))+double(0.5)), double(1))
    var $dZ#11941 : double = fmod(((($12034[int32(2)]-double(0.70456500505464))/double(0.0011352949453552))+double(0.5)), double(1))
    var $oneMinusdX#11942 : double = (double(1)-$dX#11939)
    var $oneMinusdY#11943 : double = (double(1)-$dY#11940)
    var $oneMinusdZ#11944 : double = (double(1)-$dZ#11941)
    var $weight00#11945 : double = (($12035*$oneMinusdX#11942)+($12036*$dX#11939))
    var $weight10#11946 : double = (($12037*$oneMinusdX#11942)+($12038*$dX#11939))
    var $weight01#11947 : double = (($12039*$oneMinusdX#11942)+($12040*$dX#11939))
    var $weight11#11948 : double = (($12041*$oneMinusdX#11942)+($12042*$dX#11939))
    var $weight0#11949 : double = (($weight00#11945*$oneMinusdY#11943)+($weight10#11946*$dY#11940))
    var $weight1#11950 : double = (($weight01#11947*$oneMinusdY#11943)+($weight11#11948*$dY#11940))
    $12043 = (($weight0#11949*$oneMinusdZ#11944)+($weight1#11950*$dZ#11941))
  end
  var $dummy#12044 : int32 = 0
  return $12043
end
task Particles_AddFlowCoupling($dom#11383 : region#2852(ispace#2852(int1d), particles_columns), $particles#11385 : region#2853(ispace#2853(int1d), particles_columns), $Fluid#11388 : region#2854(ispace#2854(int3d), Fluid_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#11388.centerCoordinates), reads($Fluid#11388.temperature), reads($Fluid#11388.velocity), reads($particles#11385.cell), reads($particles#11385.deltaTemperatureTerm), writes($particles#11385.deltaTemperatureTerm), reads($particles#11385.deltaVelocityOverRelaxationTime), writes($particles#11385.deltaVelocityOverRelaxationTime), reads($particles#11385.density), reads($particles#11385.diameter), reads($particles#11385.particle_temperature), reads($particles#11385.particle_velocity), reads($particles#11385.position), reads($particles#11385.position_t), writes($particles#11385.position_t), reads($particles#11385.temperature_t), writes($particles#11385.temperature_t), reads($particles#11385.velocity_t), writes($particles#11385.velocity_t), reads($particles#11385.velocity_t), writes($particles#11385.velocity_t), reads($particles#11385.__valid), $dom#11383 <= $particles#11385
do
  for $p#12157 : int1d(particles_columns, $dom#11383) in $dom#11383 do
    if $particles#11385[$p#12157].__valid then
      var $12218 : double[3]
      do
        var $12216 : int3d = $particles#11385[$p#12157].cell
        var $12217 : double[3] = $particles#11385[$p#12157].position
        var $velocity000#11645 : double[3] = double[3](array(double(0), double(0), double(0)))
        var $velocity100#11646 : double[3] = double[3](array(double(0), double(0), double(0)))
        var $velocity010#11647 : double[3] = double[3](array(double(0), double(0), double(0)))
        var $velocity110#11648 : double[3] = double[3](array(double(0), double(0), double(0)))
        var $velocity001#11649 : double[3] = double[3](array(double(0), double(0), double(0)))
        var $velocity101#11650 : double[3] = double[3](array(double(0), double(0), double(0)))
        var $velocity011#11651 : double[3] = double[3](array(double(0), double(0), double(0)))
        var $velocity111#11652 : double[3] = double[3](array(double(0), double(0), double(0)))
        var $velocity0#11653 : double[3] = $Fluid#11388[$12216].velocity
        if ($12217[int32(0)]>$Fluid#11388[$12216].centerCoordinates[int32(0)]) then
          var $velocityb#11654 : double[3] = $Fluid#11388[(($12216+{1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
          if ($12217[int32(1)]>$Fluid#11388[$12216].centerCoordinates[int32(1)]) then
            var $velocityaa#11655 : double[3] = $Fluid#11388[(($12216+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
            var $velocityab#11656 : double[3] = $Fluid#11388[(($12216+{1, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
            if ($12217[int32(2)]>$Fluid#11388[$12216].centerCoordinates[int32(2)]) then
              $velocity000#11645 = $velocity0#11653
              $velocity100#11646 = $velocityb#11654
              $velocity010#11647 = $velocityaa#11655
              $velocity110#11648 = $velocityab#11656
              $velocity001#11649 = $Fluid#11388[(($12216+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity101#11650 = $Fluid#11388[(($12216+{1, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity011#11651 = $Fluid#11388[(($12216+{0, 1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity111#11652 = $Fluid#11388[(($12216+{1, 1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
            else
              $velocity000#11645 = $Fluid#11388[(($12216+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity100#11646 = $Fluid#11388[(($12216+{1, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity010#11647 = $Fluid#11388[(($12216+{0, 1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity110#11648 = $Fluid#11388[(($12216+{1, 1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity001#11649 = $velocity0#11653
              $velocity101#11650 = $velocityb#11654
              $velocity011#11651 = $velocityaa#11655
              $velocity111#11652 = $velocityab#11656
            end
          else
            var $velocityaa#11657 : double[3] = $Fluid#11388[(($12216+{0, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
            var $velocityab#11658 : double[3] = $Fluid#11388[(($12216+{1, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
            if ($12217[int32(2)]>$Fluid#11388[$12216].centerCoordinates[int32(2)]) then
              $velocity000#11645 = $velocityaa#11657
              $velocity100#11646 = $velocityab#11658
              $velocity010#11647 = $velocity0#11653
              $velocity110#11648 = $velocityb#11654
              $velocity001#11649 = $Fluid#11388[(($12216+{0, -1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity101#11650 = $Fluid#11388[(($12216+{1, -1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity011#11651 = $Fluid#11388[(($12216+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity111#11652 = $Fluid#11388[(($12216+{1, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
            else
              $velocity000#11645 = $Fluid#11388[(($12216+{0, -1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity100#11646 = $Fluid#11388[(($12216+{1, -1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity010#11647 = $Fluid#11388[(($12216+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity110#11648 = $Fluid#11388[(($12216+{1, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity001#11649 = $velocityaa#11657
              $velocity101#11650 = $velocityab#11658
              $velocity011#11651 = $velocity0#11653
              $velocity111#11652 = $velocityb#11654
            end
          end
        else
          var $velocitya#11659 : double[3] = $Fluid#11388[(($12216+{-1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
          if ($12217[int32(1)]>$Fluid#11388[$12216].centerCoordinates[int32(1)]) then
            var $velocityaa#11660 : double[3] = $Fluid#11388[(($12216+{-1, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
            var $velocityab#11661 : double[3] = $Fluid#11388[(($12216+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
            if ($12217[int32(2)]>$Fluid#11388[$12216].centerCoordinates[int32(2)]) then
              $velocity000#11645 = $velocitya#11659
              $velocity100#11646 = $velocity0#11653
              $velocity010#11647 = $velocityaa#11660
              $velocity110#11648 = $velocityab#11661
              $velocity001#11649 = $Fluid#11388[(($12216+{-1, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity101#11650 = $Fluid#11388[(($12216+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity011#11651 = $Fluid#11388[(($12216+{-1, 1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity111#11652 = $Fluid#11388[(($12216+{0, 1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
            else
              $velocity000#11645 = $Fluid#11388[(($12216+{-1, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity100#11646 = $Fluid#11388[(($12216+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity010#11647 = $Fluid#11388[(($12216+{-1, 1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity110#11648 = $Fluid#11388[(($12216+{0, 1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity001#11649 = $velocitya#11659
              $velocity101#11650 = $velocity0#11653
              $velocity011#11651 = $velocityaa#11660
              $velocity111#11652 = $velocityab#11661
            end
          else
            var $velocityaa#11662 : double[3] = $Fluid#11388[(($12216+{-1, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
            var $velocityab#11663 : double[3] = $Fluid#11388[(($12216+{0, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
            if ($12217[int32(2)]>$Fluid#11388[$12216].centerCoordinates[int32(2)]) then
              $velocity000#11645 = $velocityaa#11662
              $velocity100#11646 = $velocityab#11663
              $velocity010#11647 = $velocitya#11659
              $velocity110#11648 = $velocity0#11653
              $velocity001#11649 = $Fluid#11388[(($12216+{-1, -1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity101#11650 = $Fluid#11388[(($12216+{0, -1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity011#11651 = $Fluid#11388[(($12216+{-1, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity111#11652 = $Fluid#11388[(($12216+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
            else
              $velocity000#11645 = $Fluid#11388[(($12216+{-1, -1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity100#11646 = $Fluid#11388[(($12216+{0, -1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity010#11647 = $Fluid#11388[(($12216+{-1, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity110#11648 = $Fluid#11388[(($12216+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].velocity
              $velocity001#11649 = $velocityaa#11662
              $velocity101#11650 = $velocityab#11663
              $velocity011#11651 = $velocitya#11659
              $velocity111#11652 = $velocity0#11653
            end
          end
        end
        var $11740 : double[3]
        do
          var $11731 : double[3] = $12217
          var $11732 : double[3] = $velocity000#11645
          var $11733 : double[3] = $velocity100#11646
          var $11734 : double[3] = $velocity010#11647
          var $11735 : double[3] = $velocity110#11648
          var $11736 : double[3] = $velocity001#11649
          var $11737 : double[3] = $velocity101#11650
          var $11738 : double[3] = $velocity011#11651
          var $11739 : double[3] = $velocity111#11652
          var $dX#11532 : double = fmod(((($11731[int32(0)]-double(0.6987052))/double(0.0015841600529101))+double(0.5)), double(1))
          var $dY#11533 : double = fmod(((($11731[int32(1)]-double(0.99700342743682))/double(0.0010188725631769))+double(0.5)), double(1))
          var $dZ#11534 : double = fmod(((($11731[int32(2)]-double(0.70456500505464))/double(0.0011352949453552))+double(0.5)), double(1))
          var $oneMinusdX#11535 : double = (double(1)-$dX#11532)
          var $oneMinusdY#11536 : double = (double(1)-$dY#11533)
          var $oneMinusdZ#11537 : double = (double(1)-$dZ#11534)
          var $weight00#11538 : double[3] = vv_add_double_3(vs_mul_double_3($11732, $oneMinusdX#11535), vs_mul_double_3($11733, $dX#11532))
          var $weight10#11539 : double[3] = vv_add_double_3(vs_mul_double_3($11734, $oneMinusdX#11535), vs_mul_double_3($11735, $dX#11532))
          var $weight01#11540 : double[3] = vv_add_double_3(vs_mul_double_3($11736, $oneMinusdX#11535), vs_mul_double_3($11737, $dX#11532))
          var $weight11#11541 : double[3] = vv_add_double_3(vs_mul_double_3($11738, $oneMinusdX#11535), vs_mul_double_3($11739, $dX#11532))
          var $weight0#11542 : double[3] = vv_add_double_3(vs_mul_double_3($weight00#11538, $oneMinusdY#11536), vs_mul_double_3($weight10#11539, $dY#11533))
          var $weight1#11543 : double[3] = vv_add_double_3(vs_mul_double_3($weight01#11540, $oneMinusdY#11536), vs_mul_double_3($weight11#11541, $dY#11533))
          $11740 = vv_add_double_3(vs_mul_double_3($weight0#11542, $oneMinusdZ#11537), vs_mul_double_3($weight1#11543, $dZ#11534))
        end
        var $dummy#11741 : int32 = 0
        $12218 = $11740
      end
      var $dummy#12219 : int32 = 0
      var $flowVelocity#12158 : double[3] = $12218
      var $12222 : double
      do
        var $12220 : int3d = $particles#11385[$p#12157].cell
        var $12221 : double[3] = $particles#11385[$p#12157].position
        var $temp000#11996 : double = double(double(0))
        var $temp100#11997 : double = double(double(0))
        var $temp010#11998 : double = double(double(0))
        var $temp110#11999 : double = double(double(0))
        var $temp001#12000 : double = double(double(0))
        var $temp101#12001 : double = double(double(0))
        var $temp011#12002 : double = double(double(0))
        var $temp111#12003 : double = double(double(0))
        var $temp0#12004 : double = $Fluid#11388[$12220].temperature
        if ($12221[int32(0)]>$Fluid#11388[$12220].centerCoordinates[int32(0)]) then
          var $tempb#12005 : double = $Fluid#11388[(($12220+{1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
          if ($12221[int32(1)]>$Fluid#11388[$12220].centerCoordinates[int32(1)]) then
            var $tempaa#12006 : double = $Fluid#11388[(($12220+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
            var $tempab#12007 : double = $Fluid#11388[(($12220+{1, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
            if ($12221[int32(2)]>$Fluid#11388[$12220].centerCoordinates[int32(2)]) then
              $temp000#11996 = $temp0#12004
              $temp100#11997 = $tempb#12005
              $temp010#11998 = $tempaa#12006
              $temp110#11999 = $tempab#12007
              $temp001#12000 = $Fluid#11388[(($12220+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp101#12001 = $Fluid#11388[(($12220+{1, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp011#12002 = $Fluid#11388[(($12220+{0, 1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp111#12003 = $Fluid#11388[(($12220+{1, 1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
            else
              $temp000#11996 = $Fluid#11388[(($12220+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp100#11997 = $Fluid#11388[(($12220+{1, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp010#11998 = $Fluid#11388[(($12220+{0, 1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp110#11999 = $Fluid#11388[(($12220+{1, 1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp001#12000 = $temp0#12004
              $temp101#12001 = $tempb#12005
              $temp011#12002 = $tempaa#12006
              $temp111#12003 = $tempab#12007
            end
          else
            var $tempaa#12008 : double = $Fluid#11388[(($12220+{0, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
            var $tempab#12009 : double = $Fluid#11388[(($12220+{1, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
            if ($12221[int32(2)]>$Fluid#11388[$12220].centerCoordinates[int32(2)]) then
              $temp000#11996 = $tempaa#12008
              $temp100#11997 = $tempab#12009
              $temp010#11998 = $temp0#12004
              $temp110#11999 = $tempb#12005
              $temp001#12000 = $Fluid#11388[(($12220+{0, -1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp101#12001 = $Fluid#11388[(($12220+{1, -1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp011#12002 = $Fluid#11388[(($12220+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp111#12003 = $Fluid#11388[(($12220+{1, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
            else
              $temp000#11996 = $Fluid#11388[(($12220+{0, -1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp100#11997 = $Fluid#11388[(($12220+{1, -1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp010#11998 = $Fluid#11388[(($12220+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp110#11999 = $Fluid#11388[(($12220+{1, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp001#12000 = $tempaa#12008
              $temp101#12001 = $tempab#12009
              $temp011#12002 = $temp0#12004
              $temp111#12003 = $tempb#12005
            end
          end
        else
          var $tempa#12010 : double = $Fluid#11388[(($12220+{-1, 0, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
          if ($12221[int32(1)]>$Fluid#11388[$12220].centerCoordinates[int32(1)]) then
            var $tempaa#12011 : double = $Fluid#11388[(($12220+{-1, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
            var $tempab#12012 : double = $Fluid#11388[(($12220+{0, 1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
            if ($12221[int32(2)]>$Fluid#11388[$12220].centerCoordinates[int32(2)]) then
              $temp000#11996 = $tempa#12010
              $temp100#11997 = $temp0#12004
              $temp010#11998 = $tempaa#12011
              $temp110#11999 = $tempab#12012
              $temp001#12000 = $Fluid#11388[(($12220+{-1, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp101#12001 = $Fluid#11388[(($12220+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp011#12002 = $Fluid#11388[(($12220+{-1, 1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp111#12003 = $Fluid#11388[(($12220+{0, 1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
            else
              $temp000#11996 = $Fluid#11388[(($12220+{-1, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp100#11997 = $Fluid#11388[(($12220+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp010#11998 = $Fluid#11388[(($12220+{-1, 1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp110#11999 = $Fluid#11388[(($12220+{0, 1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp001#12000 = $tempa#12010
              $temp101#12001 = $temp0#12004
              $temp011#12002 = $tempaa#12011
              $temp111#12003 = $tempab#12012
            end
          else
            var $tempaa#12013 : double = $Fluid#11388[(($12220+{-1, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
            var $tempab#12014 : double = $Fluid#11388[(($12220+{0, -1, 0})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
            if ($12221[int32(2)]>$Fluid#11388[$12220].centerCoordinates[int32(2)]) then
              $temp000#11996 = $tempaa#12013
              $temp100#11997 = $tempab#12014
              $temp010#11998 = $tempa#12010
              $temp110#11999 = $temp0#12004
              $temp001#12000 = $Fluid#11388[(($12220+{-1, -1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp101#12001 = $Fluid#11388[(($12220+{0, -1, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp011#12002 = $Fluid#11388[(($12220+{-1, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp111#12003 = $Fluid#11388[(($12220+{0, 0, 1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
            else
              $temp000#11996 = $Fluid#11388[(($12220+{-1, -1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp100#11997 = $Fluid#11388[(($12220+{0, -1, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp010#11998 = $Fluid#11388[(($12220+{-1, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp110#11999 = $Fluid#11388[(($12220+{0, 0, -1})%rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 832, z = 733})}))].temperature
              $temp001#12000 = $tempaa#12013
              $temp101#12001 = $tempab#12014
              $temp011#12002 = $tempa#12010
              $temp111#12003 = $temp0#12004
            end
          end
        end
        var $12043 : double
        do
          var $12034 : double[3] = $12221
          var $12035 : double = $temp000#11996
          var $12036 : double = $temp100#11997
          var $12037 : double = $temp010#11998
          var $12038 : double = $temp110#11999
          var $12039 : double = $temp001#12000
          var $12040 : double = $temp101#12001
          var $12041 : double = $temp011#12002
          var $12042 : double = $temp111#12003
          var $dX#11939 : double = fmod(((($12034[int32(0)]-double(0.6987052))/double(0.0015841600529101))+double(0.5)), double(1))
          var $dY#11940 : double = fmod(((($12034[int32(1)]-double(0.99700342743682))/double(0.0010188725631769))+double(0.5)), double(1))
          var $dZ#11941 : double = fmod(((($12034[int32(2)]-double(0.70456500505464))/double(0.0011352949453552))+double(0.5)), double(1))
          var $oneMinusdX#11942 : double = (double(1)-$dX#11939)
          var $oneMinusdY#11943 : double = (double(1)-$dY#11940)
          var $oneMinusdZ#11944 : double = (double(1)-$dZ#11941)
          var $weight00#11945 : double = (($12035*$oneMinusdX#11942)+($12036*$dX#11939))
          var $weight10#11946 : double = (($12037*$oneMinusdX#11942)+($12038*$dX#11939))
          var $weight01#11947 : double = (($12039*$oneMinusdX#11942)+($12040*$dX#11939))
          var $weight11#11948 : double = (($12041*$oneMinusdX#11942)+($12042*$dX#11939))
          var $weight0#11949 : double = (($weight00#11945*$oneMinusdY#11943)+($weight10#11946*$dY#11940))
          var $weight1#11950 : double = (($weight01#11947*$oneMinusdY#11943)+($weight11#11948*$dY#11940))
          $12043 = (($weight0#11949*$oneMinusdZ#11944)+($weight1#11950*$dZ#11941))
        end
        var $dummy#12044 : int32 = 0
        $12222 = $12043
      end
      var $dummy#12223 : int32 = 0
      var $flowTemperature#12159 : double = $12222
      var $12225 : double
      do
        var $12224 : double = $flowTemperature#12159
        var $viscosity#5391 : double = double(double(0))
        if (int32(3003)==int32(3001)) then
          $viscosity#5391 = double(0.2454194)
        else
          if (int32(3003)==int32(3002)) then
            $viscosity#5391 = (double(0.9683056)*pow(($12224/double(0.9913615)), double(0.75)))
          else
            if (int32(3003)==int32(3003)) then
              $viscosity#5391 = ((double(0.5815515)*pow(($12224/double(0.6935851)), (double(3)/double(2))))*((double(0.6935851)+double(0.2772065))/($12224+double(0.2772065))))
            else
              std.assert(false, "(Liszt assertion)")
            end
          end
        end
        $12225 = $viscosity#5391
      end
      var $dummy#12226 : int32 = 0
      var $flowDynamicViscosity#12160 : double = $12225
      if (int32(6002)==int32(6001)) then
      else
        if (int32(6002)==int32(6002)) then
          var $tmp#12161 : double[3] = $particles#11385[$p#12157].particle_velocity
          var $v#12162 : double[3] = $particles#11385[$p#12157].position_t
          $v#12162[0] += $tmp#12161[0]
          $v#12162[1] += $tmp#12161[1]
          $v#12162[2] += $tmp#12161[2]
          $particles#11385[$p#12157].position_t = $v#12162
        else
          std.assert(false, "(Liszt assertion)")
        end
      end
      var $particleReynoldsNumber#12163 : double = double(0)
      var $relaxationTime#12164 : double = ((($particles#11385[$p#12157].density*pow($particles#11385[$p#12157].diameter, double(int32(2))))/(double(18)*$flowDynamicViscosity#12160))/(double(1)+(double(0.15)*pow($particleReynoldsNumber#12163, double(0.687)))))
      $particles#11385[$p#12157].deltaVelocityOverRelaxationTime = vs_div_double_3(vv_sub_double_3($flowVelocity#12158, $particles#11385[$p#12157].particle_velocity), $relaxationTime#12164)
      $particles#11385[$p#12157].deltaTemperatureTerm = (((double(3.1415926535898)*pow($particles#11385[$p#12157].diameter, double(int32(2))))*double(0.8720309))*($flowTemperature#12159-$particles#11385[$p#12157].particle_temperature))
      if (int32(6002)==int32(6001)) then
        $particles#11385[$p#12157].velocity_t = array(double(0), double(0), double(0))
      else
        if (int32(6002)==int32(6002)) then
          var $tmp#12165 : double[3] = $particles#11385[$p#12157].deltaVelocityOverRelaxationTime
          var $v#12166 : double[3] = $particles#11385[$p#12157].velocity_t
          $v#12166[0] += $tmp#12165[0]
          $v#12166[1] += $tmp#12165[1]
          $v#12166[2] += $tmp#12165[2]
          $particles#11385[$p#12157].velocity_t = $v#12166
        else
          std.assert(false, "(Liszt assertion)")
        end
      end
      $particles#11385[$p#12157].temperature_t += ($particles#11385[$p#12157].deltaTemperatureTerm/((((double(3.1415926535898)*pow($particles#11385[$p#12157].diameter, double(int32(3))))/double(6))*$particles#11385[$p#12157].density)*double(0.3767944)))
    else
    end
  end
end
task Particles_AddBodyForces($dom#13198 : region#3527(ispace#3527(int1d), particles_columns), $particles#13200 : region#3528(ispace#3528(int1d), particles_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($particles#13200.velocity_t), writes($particles#13200.velocity_t), reads($particles#13200.__valid), $dom#13198 <= $particles#13200
do
  for $p#13215 : int1d(particles_columns, $dom#13198) in $dom#13198 do
    if $particles#13200[$p#13215].__valid then
      var $tmp#13216 : double[3] = array(double(0.5899043), double(0.3012763), double(0.6355317))
      var $v#13217 : double[3] = $particles#13200[$p#13215].velocity_t
      $v#13217[0] += $tmp#13216[0]
      $v#13217[1] += $tmp#13216[1]
      $v#13217[2] += $tmp#13216[2]
      $particles#13200[$p#13215].velocity_t = $v#13217
    else
    end
  end
end
task Flow_AddParticlesCoupling($dom#13253 : region#3542(ispace#3542(int1d), particles_columns), $particles#13255 : region#3543(ispace#3543(int1d), particles_columns), $Fluid#13258 : region#3544(ispace#3544(int3d), Fluid_columns)) : double
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#13258.rhoEnergy_t), writes($Fluid#13258.rhoEnergy_t), reads($Fluid#13258.rhoVelocity_t), writes($Fluid#13258.rhoVelocity_t), reads($particles#13255.cell), reads($particles#13255.deltaTemperatureTerm), reads($particles#13255.deltaVelocityOverRelaxationTime), reads($particles#13255.density), reads($particles#13255.diameter), reads($particles#13255.__valid), $dom#13253 <= $particles#13255
do
  var $acc#13279 : double = double(0)
  for $p#13280 : int1d(particles_columns, $dom#13253) in $dom#13253 do
    if $particles#13255[$p#13280].__valid then
      var $cellVolume#13281 : double = ((double(0.0015841600529101)*double(0.0010188725631769))*double(0.0011352949453552))
      var $tmp#13282 : double[3] = vs_div_double_3(vs_mul_double_3($particles#13255[$p#13280].deltaVelocityOverRelaxationTime, (-(((double(3.1415926535898)*pow($particles#13255[$p#13280].diameter, double(int32(3))))/double(6))*$particles#13255[$p#13280].density))), $cellVolume#13281)
      var $v#13283 : double[3] = $Fluid#13258[$particles#13255[$p#13280].cell].rhoVelocity_t
      $v#13283[0] += $tmp#13282[0]
      $v#13283[1] += $tmp#13282[1]
      $v#13283[2] += $tmp#13282[2]
      $Fluid#13258[$particles#13255[$p#13280].cell].rhoVelocity_t = $v#13283
      $Fluid#13258[$particles#13255[$p#13280].cell].rhoEnergy_t += ((-$particles#13255[$p#13280].deltaTemperatureTerm)/$cellVolume#13281)
      if false then
        $acc#13279 += ($particles#13255[$p#13280].deltaTemperatureTerm/$cellVolume#13281)
      else
      end
    else
    end
  end
  return $acc#13279
end
task Flow_UpdateVars($dom#13344 : region#3568(ispace#3568(int3d), Fluid_columns), $Fluid#13346 : region#3569(ispace#3569(int3d), Fluid_columns), $TimeIntegrator_deltaTime#13348 : double, $TimeIntegrator_stage#13349 : int32)
-- leaf (false), inner (false), idempotent (false)
where
  reads($Fluid#13346.rho), writes($Fluid#13346.rho), reads($Fluid#13346.rhoEnergy), writes($Fluid#13346.rhoEnergy), reads($Fluid#13346.rhoEnergy_new), reads($Fluid#13346.rhoEnergy_new), writes($Fluid#13346.rhoEnergy_new), reads($Fluid#13346.rhoEnergy_old), reads($Fluid#13346.rhoEnergy_t), reads($Fluid#13346.rhoVelocity), writes($Fluid#13346.rhoVelocity), reads($Fluid#13346.rhoVelocity_new), reads($Fluid#13346.rhoVelocity_new), writes($Fluid#13346.rhoVelocity_new), reads($Fluid#13346.rhoVelocity_old), reads($Fluid#13346.rhoVelocity_t), reads($Fluid#13346.rho_new), reads($Fluid#13346.rho_new), writes($Fluid#13346.rho_new), reads($Fluid#13346.rho_old), reads($Fluid#13346.rho_t), $dom#13344 <= $Fluid#13346
do
  for $c#13403 : int3d(Fluid_columns, $dom#13344) in $dom#13344 do
    var $deltaTime#13404 : double = $TimeIntegrator_deltaTime#13348
    if ($TimeIntegrator_stage#13349==int32(1)) then
      $Fluid#13346[$c#13403].rho_new += (((double(1)/double(6))*$deltaTime#13404)*$Fluid#13346[$c#13403].rho_t)
      $Fluid#13346[$c#13403].rho = ($Fluid#13346[$c#13403].rho_old+((double(0.5)*$deltaTime#13404)*$Fluid#13346[$c#13403].rho_t))
      var $tmp#13405 : double[3] = vs_mul_double_3($Fluid#13346[$c#13403].rhoVelocity_t, ((double(1)/double(6))*$deltaTime#13404))
      var $v#13406 : double[3] = $Fluid#13346[$c#13403].rhoVelocity_new
      $v#13406[0] += $tmp#13405[0]
      $v#13406[1] += $tmp#13405[1]
      $v#13406[2] += $tmp#13405[2]
      $Fluid#13346[$c#13403].rhoVelocity_new = $v#13406
      $Fluid#13346[$c#13403].rhoVelocity = vv_add_double_3($Fluid#13346[$c#13403].rhoVelocity_old, vs_mul_double_3($Fluid#13346[$c#13403].rhoVelocity_t, (double(0.5)*$deltaTime#13404)))
      $Fluid#13346[$c#13403].rhoEnergy_new += (((double(1)/double(6))*$deltaTime#13404)*$Fluid#13346[$c#13403].rhoEnergy_t)
      $Fluid#13346[$c#13403].rhoEnergy = ($Fluid#13346[$c#13403].rhoEnergy_old+((double(0.5)*$deltaTime#13404)*$Fluid#13346[$c#13403].rhoEnergy_t))
    else
      if ($TimeIntegrator_stage#13349==int32(2)) then
        $Fluid#13346[$c#13403].rho_new += (((double(1)/double(3))*$deltaTime#13404)*$Fluid#13346[$c#13403].rho_t)
        $Fluid#13346[$c#13403].rho = ($Fluid#13346[$c#13403].rho_old+((double(0.5)*$deltaTime#13404)*$Fluid#13346[$c#13403].rho_t))
        var $tmp#13407 : double[3] = vs_mul_double_3($Fluid#13346[$c#13403].rhoVelocity_t, ((double(1)/double(3))*$deltaTime#13404))
        var $v#13408 : double[3] = $Fluid#13346[$c#13403].rhoVelocity_new
        $v#13408[0] += $tmp#13407[0]
        $v#13408[1] += $tmp#13407[1]
        $v#13408[2] += $tmp#13407[2]
        $Fluid#13346[$c#13403].rhoVelocity_new = $v#13408
        $Fluid#13346[$c#13403].rhoVelocity = vv_add_double_3($Fluid#13346[$c#13403].rhoVelocity_old, vs_mul_double_3($Fluid#13346[$c#13403].rhoVelocity_t, (double(0.5)*$deltaTime#13404)))
        $Fluid#13346[$c#13403].rhoEnergy_new += (((double(1)/double(3))*$deltaTime#13404)*$Fluid#13346[$c#13403].rhoEnergy_t)
        $Fluid#13346[$c#13403].rhoEnergy = ($Fluid#13346[$c#13403].rhoEnergy_old+((double(0.5)*$deltaTime#13404)*$Fluid#13346[$c#13403].rhoEnergy_t))
      else
        if ($TimeIntegrator_stage#13349==int32(3)) then
          $Fluid#13346[$c#13403].rho_new += (((double(1)/double(3))*$deltaTime#13404)*$Fluid#13346[$c#13403].rho_t)
          $Fluid#13346[$c#13403].rho = ($Fluid#13346[$c#13403].rho_old+((double(1)*$deltaTime#13404)*$Fluid#13346[$c#13403].rho_t))
          var $tmp#13409 : double[3] = vs_mul_double_3($Fluid#13346[$c#13403].rhoVelocity_t, ((double(1)/double(3))*$deltaTime#13404))
          var $v#13410 : double[3] = $Fluid#13346[$c#13403].rhoVelocity_new
          $v#13410[0] += $tmp#13409[0]
          $v#13410[1] += $tmp#13409[1]
          $v#13410[2] += $tmp#13409[2]
          $Fluid#13346[$c#13403].rhoVelocity_new = $v#13410
          $Fluid#13346[$c#13403].rhoVelocity = vv_add_double_3($Fluid#13346[$c#13403].rhoVelocity_old, vs_mul_double_3($Fluid#13346[$c#13403].rhoVelocity_t, (double(1)*$deltaTime#13404)))
          $Fluid#13346[$c#13403].rhoEnergy_new += (((double(1)/double(3))*$deltaTime#13404)*$Fluid#13346[$c#13403].rhoEnergy_t)
          $Fluid#13346[$c#13403].rhoEnergy = ($Fluid#13346[$c#13403].rhoEnergy_old+((double(1)*$deltaTime#13404)*$Fluid#13346[$c#13403].rhoEnergy_t))
        else
          $Fluid#13346[$c#13403].rho = ($Fluid#13346[$c#13403].rho_new+(((double(1)/double(6))*$deltaTime#13404)*$Fluid#13346[$c#13403].rho_t))
          $Fluid#13346[$c#13403].rhoVelocity = vv_add_double_3($Fluid#13346[$c#13403].rhoVelocity_new, vs_mul_double_3($Fluid#13346[$c#13403].rhoVelocity_t, ((double(1)/double(6))*$deltaTime#13404)))
          $Fluid#13346[$c#13403].rhoEnergy = ($Fluid#13346[$c#13403].rhoEnergy_new+(((double(1)/double(6))*$deltaTime#13404)*$Fluid#13346[$c#13403].rhoEnergy_t))
        end
      end
    end
  end
end
task Particles_UpdateVars($dom#13503 : region#3597(ispace#3597(int1d), particles_columns), $particles#13505 : region#3598(ispace#3598(int1d), particles_columns), $TimeIntegrator_deltaTime#13507 : double, $TimeIntegrator_stage#13508 : int32)
-- leaf (false), inner (false), idempotent (false)
where
  reads($particles#13505.particle_temperature), writes($particles#13505.particle_temperature), reads($particles#13505.particle_velocity), writes($particles#13505.particle_velocity), reads($particles#13505.position), writes($particles#13505.position), reads($particles#13505.position_new), reads($particles#13505.position_new), writes($particles#13505.position_new), reads($particles#13505.position_old), reads($particles#13505.position_t), reads($particles#13505.temperature_new), reads($particles#13505.temperature_new), writes($particles#13505.temperature_new), reads($particles#13505.temperature_old), reads($particles#13505.temperature_t), reads($particles#13505.velocity_new), reads($particles#13505.velocity_new), writes($particles#13505.velocity_new), reads($particles#13505.velocity_old), reads($particles#13505.velocity_t), reads($particles#13505.__valid), $dom#13503 <= $particles#13505
do
  for $p#13623 : int1d(particles_columns, $dom#13503) in $dom#13503 do
    if $particles#13505[$p#13623].__valid then
      var $deltaTime#13624 : double = $TimeIntegrator_deltaTime#13507
      if ($TimeIntegrator_stage#13508==int32(1)) then
        var $tmp#13625 : double[3] = vs_mul_double_3($particles#13505[$p#13623].position_t, ((double(1)/double(6))*$deltaTime#13624))
        var $v#13626 : double[3] = $particles#13505[$p#13623].position_new
        $v#13626[0] += $tmp#13625[0]
        $v#13626[1] += $tmp#13625[1]
        $v#13626[2] += $tmp#13625[2]
        $particles#13505[$p#13623].position_new = $v#13626
        $particles#13505[$p#13623].position = vv_add_double_3($particles#13505[$p#13623].position_old, vs_mul_double_3($particles#13505[$p#13623].position_t, (double(0.5)*$deltaTime#13624)))
        var $tmp#13627 : double[3] = vs_mul_double_3($particles#13505[$p#13623].velocity_t, ((double(1)/double(6))*$deltaTime#13624))
        var $v#13628 : double[3] = $particles#13505[$p#13623].velocity_new
        $v#13628[0] += $tmp#13627[0]
        $v#13628[1] += $tmp#13627[1]
        $v#13628[2] += $tmp#13627[2]
        $particles#13505[$p#13623].velocity_new = $v#13628
        $particles#13505[$p#13623].particle_velocity = vv_add_double_3($particles#13505[$p#13623].velocity_old, vs_mul_double_3($particles#13505[$p#13623].velocity_t, (double(0.5)*$deltaTime#13624)))
        $particles#13505[$p#13623].temperature_new += (((double(1)/double(6))*$deltaTime#13624)*$particles#13505[$p#13623].temperature_t)
        $particles#13505[$p#13623].particle_temperature = ($particles#13505[$p#13623].temperature_old+((double(0.5)*$deltaTime#13624)*$particles#13505[$p#13623].temperature_t))
      else
        if ($TimeIntegrator_stage#13508==int32(2)) then
          var $tmp#13629 : double[3] = vs_mul_double_3($particles#13505[$p#13623].position_t, ((double(1)/double(3))*$deltaTime#13624))
          var $v#13630 : double[3] = $particles#13505[$p#13623].position_new
          $v#13630[0] += $tmp#13629[0]
          $v#13630[1] += $tmp#13629[1]
          $v#13630[2] += $tmp#13629[2]
          $particles#13505[$p#13623].position_new = $v#13630
          $particles#13505[$p#13623].position = vv_add_double_3($particles#13505[$p#13623].position_old, vs_mul_double_3($particles#13505[$p#13623].position_t, (double(0.5)*$deltaTime#13624)))
          var $tmp#13631 : double[3] = vs_mul_double_3($particles#13505[$p#13623].velocity_t, ((double(1)/double(3))*$deltaTime#13624))
          var $v#13632 : double[3] = $particles#13505[$p#13623].velocity_new
          $v#13632[0] += $tmp#13631[0]
          $v#13632[1] += $tmp#13631[1]
          $v#13632[2] += $tmp#13631[2]
          $particles#13505[$p#13623].velocity_new = $v#13632
          $particles#13505[$p#13623].particle_velocity = vv_add_double_3($particles#13505[$p#13623].velocity_old, vs_mul_double_3($particles#13505[$p#13623].velocity_t, (double(0.5)*$deltaTime#13624)))
          $particles#13505[$p#13623].temperature_new += (((double(1)/double(3))*$deltaTime#13624)*$particles#13505[$p#13623].temperature_t)
          $particles#13505[$p#13623].particle_temperature = ($particles#13505[$p#13623].temperature_old+((double(0.5)*$deltaTime#13624)*$particles#13505[$p#13623].temperature_t))
        else
          if ($TimeIntegrator_stage#13508==int32(3)) then
            var $tmp#13633 : double[3] = vs_mul_double_3($particles#13505[$p#13623].position_t, ((double(1)/double(3))*$deltaTime#13624))
            var $v#13634 : double[3] = $particles#13505[$p#13623].position_new
            $v#13634[0] += $tmp#13633[0]
            $v#13634[1] += $tmp#13633[1]
            $v#13634[2] += $tmp#13633[2]
            $particles#13505[$p#13623].position_new = $v#13634
            $particles#13505[$p#13623].position = vv_add_double_3($particles#13505[$p#13623].position_old, vs_mul_double_3($particles#13505[$p#13623].position_t, (double(1)*$deltaTime#13624)))
            var $tmp#13635 : double[3] = vs_mul_double_3($particles#13505[$p#13623].velocity_t, ((double(1)/double(3))*$deltaTime#13624))
            var $v#13636 : double[3] = $particles#13505[$p#13623].velocity_new
            $v#13636[0] += $tmp#13635[0]
            $v#13636[1] += $tmp#13635[1]
            $v#13636[2] += $tmp#13635[2]
            $particles#13505[$p#13623].velocity_new = $v#13636
            $particles#13505[$p#13623].particle_velocity = vv_add_double_3($particles#13505[$p#13623].velocity_old, vs_mul_double_3($particles#13505[$p#13623].velocity_t, (double(1)*$deltaTime#13624)))
            $particles#13505[$p#13623].temperature_new += (((double(1)/double(3))*$deltaTime#13624)*$particles#13505[$p#13623].temperature_t)
            $particles#13505[$p#13623].particle_temperature = ($particles#13505[$p#13623].temperature_old+((double(1)*$deltaTime#13624)*$particles#13505[$p#13623].temperature_t))
          else
            $particles#13505[$p#13623].position = vv_add_double_3($particles#13505[$p#13623].position_new, vs_mul_double_3($particles#13505[$p#13623].position_t, ((double(1)/double(6))*$deltaTime#13624)))
            $particles#13505[$p#13623].particle_velocity = vv_add_double_3($particles#13505[$p#13623].velocity_new, vs_mul_double_3($particles#13505[$p#13623].velocity_t, ((double(1)/double(6))*$deltaTime#13624)))
            $particles#13505[$p#13623].particle_temperature = ($particles#13505[$p#13623].temperature_new+(((double(1)/double(6))*$deltaTime#13624)*$particles#13505[$p#13623].temperature_t))
          end
        end
      end
    else
    end
  end
end
task Particles_UpdateAuxiliaryStep1($dom#13800 : region#3638(ispace#3638(int1d), particles_columns), $particles#13802 : region#3639(ispace#3639(int1d), particles_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($particles#13802.particle_velocity), reads($particles#13802.position), reads($particles#13802.position_ghost), writes($particles#13802.position_ghost), reads($particles#13802.velocity_ghost), writes($particles#13802.velocity_ghost), reads($particles#13802.velocity_ghost), writes($particles#13802.velocity_ghost), reads($particles#13802.velocity_t), reads($particles#13802.velocity_t_ghost), writes($particles#13802.velocity_t_ghost), reads($particles#13802.velocity_t_ghost), writes($particles#13802.velocity_t_ghost), reads($particles#13802.__valid), $dom#13800 <= $particles#13802
do
  for $p#13925 : int1d(particles_columns, $dom#13800) in $dom#13800 do
    if $particles#13802[$p#13925].__valid then
      $particles#13802[$p#13925].position_ghost[int32(0)] = $particles#13802[$p#13925].position[int32(0)]
      $particles#13802[$p#13925].position_ghost[int32(1)] = $particles#13802[$p#13925].position[int32(1)]
      $particles#13802[$p#13925].position_ghost[int32(2)] = $particles#13802[$p#13925].position[int32(2)]
      $particles#13802[$p#13925].velocity_ghost[int32(0)] = $particles#13802[$p#13925].particle_velocity[int32(0)]
      $particles#13802[$p#13925].velocity_ghost[int32(1)] = $particles#13802[$p#13925].particle_velocity[int32(1)]
      $particles#13802[$p#13925].velocity_ghost[int32(2)] = $particles#13802[$p#13925].particle_velocity[int32(2)]
      $particles#13802[$p#13925].velocity_t_ghost[int32(0)] = $particles#13802[$p#13925].velocity_t[int32(0)]
      $particles#13802[$p#13925].velocity_t_ghost[int32(1)] = $particles#13802[$p#13925].velocity_t[int32(1)]
      $particles#13802[$p#13925].velocity_t_ghost[int32(2)] = $particles#13802[$p#13925].velocity_t[int32(2)]
      if ($particles#13802[$p#13925].position[int32(0)]<double(0.6987052)) then
        if (int32(2001)==int32(2001)) then
          $particles#13802[$p#13925].position_ghost[int32(0)] = ($particles#13802[$p#13925].position[int32(0)]+double(0.5988125))
        else
          if (int32(2001)==int32(2002)) then
            $particles#13802[$p#13925].position_ghost[int32(0)] = double(0.6987052)
            var $impulse#13926 : double = ((-(double(1)+double(0.1487899)))*$particles#13802[$p#13925].particle_velocity[int32(0)])
            if ($impulse#13926<=double(int32(0))) then
              $particles#13802[$p#13925].velocity_ghost[int32(0)] += $impulse#13926
            else
            end
            var $contact_force#13927 : double = (double(-1)*$particles#13802[$p#13925].velocity_t[int32(0)])
            if ($contact_force#13927>double(int32(0))) then
              $particles#13802[$p#13925].velocity_t_ghost[int32(0)] += $contact_force#13927
            else
            end
          else
            std.assert(false, "(Liszt assertion)")
          end
        end
      else
      end
      if ($particles#13802[$p#13925].position[int32(0)]>(double(0.6987052)+double(0.5988125))) then
        if (int32(2001)==int32(2001)) then
          $particles#13802[$p#13925].position_ghost[int32(0)] = ($particles#13802[$p#13925].position[int32(0)]-double(0.5988125))
        else
          if (int32(2001)==int32(2002)) then
            $particles#13802[$p#13925].position_ghost[int32(0)] = (double(0.6987052)+double(0.5988125))
            var $impulse#13928 : double = ((-(double(1)+double(0.1487899)))*$particles#13802[$p#13925].particle_velocity[int32(0)])
            if ($impulse#13928>=double(int32(0))) then
              $particles#13802[$p#13925].velocity_ghost[int32(0)] += $impulse#13928
            else
            end
            var $contact_force#13929 : double = (double(-1)*$particles#13802[$p#13925].velocity_t[int32(0)])
            if ($contact_force#13929<double(int32(0))) then
              $particles#13802[$p#13925].velocity_t_ghost[int32(0)] += $contact_force#13929
            else
            end
          else
            std.assert(false, "(Liszt assertion)")
          end
        end
      else
      end
      if ($particles#13802[$p#13925].position[int32(1)]<double(0.9980223)) then
        if (int32(2002)==int32(2001)) then
          $particles#13802[$p#13925].position_ghost[int32(1)] = ($particles#13802[$p#13925].position[int32(1)]+double(0.8466831))
        else
          if (int32(2002)==int32(2002)) then
            $particles#13802[$p#13925].position_ghost[int32(1)] = double(0.9980223)
            var $impulse#13930 : double = ((-(double(1)+double(0.1487899)))*$particles#13802[$p#13925].particle_velocity[int32(1)])
            if ($impulse#13930<=double(int32(0))) then
              $particles#13802[$p#13925].velocity_ghost[int32(1)] += $impulse#13930
            else
            end
            var $contact_force#13931 : double = (double(-1)*$particles#13802[$p#13925].velocity_t[int32(1)])
            if ($contact_force#13931>double(int32(0))) then
              $particles#13802[$p#13925].velocity_t_ghost[int32(1)] += $contact_force#13931
            else
            end
          else
            std.assert(false, "(Liszt assertion)")
          end
        end
      else
      end
      if ($particles#13802[$p#13925].position[int32(1)]>(double(0.9980223)+double(0.8466831))) then
        if (int32(2002)==int32(2001)) then
          $particles#13802[$p#13925].position_ghost[int32(1)] = ($particles#13802[$p#13925].position[int32(1)]-double(0.8466831))
        else
          if (int32(2002)==int32(2002)) then
            $particles#13802[$p#13925].position_ghost[int32(1)] = (double(0.9980223)+double(0.8466831))
            var $impulse#13932 : double = ((-(double(1)+double(0.1487899)))*$particles#13802[$p#13925].particle_velocity[int32(1)])
            if ($impulse#13932>=double(int32(0))) then
              $particles#13802[$p#13925].velocity_ghost[int32(1)] += $impulse#13932
            else
            end
            var $contact_force#13933 : double = (double(-1)*$particles#13802[$p#13925].velocity_t[int32(1)])
            if ($contact_force#13933<double(int32(0))) then
              $particles#13802[$p#13925].velocity_t_ghost[int32(1)] += $contact_force#13933
            else
            end
          else
            std.assert(false, "(Liszt assertion)")
          end
        end
      else
      end
      if ($particles#13802[$p#13925].position[int32(2)]<double(0.7057003)) then
        if (int32(2002)==int32(2001)) then
          $particles#13802[$p#13925].position_ghost[int32(2)] = ($particles#13802[$p#13925].position[int32(2)]+double(0.8310359))
        else
          if (int32(2002)==int32(2002)) then
            $particles#13802[$p#13925].position_ghost[int32(2)] = double(0.7057003)
            var $impulse#13934 : double = ((-(double(1)+double(0.1487899)))*$particles#13802[$p#13925].particle_velocity[int32(2)])
            if ($impulse#13934<=double(int32(0))) then
              $particles#13802[$p#13925].velocity_ghost[int32(2)] += $impulse#13934
            else
            end
            var $contact_force#13935 : double = (double(-1)*$particles#13802[$p#13925].velocity_t[int32(2)])
            if ($contact_force#13935>double(int32(0))) then
              $particles#13802[$p#13925].velocity_t_ghost[int32(2)] += $contact_force#13935
            else
            end
          else
            std.assert(false, "(Liszt assertion)")
          end
        end
      else
      end
      if ($particles#13802[$p#13925].position[int32(2)]>(double(0.7057003)+double(0.8310359))) then
        if (int32(2002)==int32(2001)) then
          $particles#13802[$p#13925].position_ghost[int32(2)] = ($particles#13802[$p#13925].position[int32(2)]-double(0.8310359))
        else
          if (int32(2002)==int32(2002)) then
            $particles#13802[$p#13925].position_ghost[int32(2)] = (double(0.7057003)+double(0.8310359))
            var $impulse#13936 : double = ((-(double(1)+double(0.1487899)))*$particles#13802[$p#13925].particle_velocity[int32(2)])
            if ($impulse#13936>=double(int32(0))) then
              $particles#13802[$p#13925].velocity_ghost[int32(2)] += $impulse#13936
            else
            end
            var $contact_force#13937 : double = (double(-1)*$particles#13802[$p#13925].velocity_t[int32(2)])
            if ($contact_force#13937<double(int32(0))) then
              $particles#13802[$p#13925].velocity_t_ghost[int32(2)] += $contact_force#13937
            else
            end
          else
            std.assert(false, "(Liszt assertion)")
          end
        end
      else
      end
    else
    end
  end
end
task Particles_UpdateAuxiliaryStep2($dom#14021 : region#3673(ispace#3673(int1d), particles_columns), $particles#14023 : region#3674(ispace#3674(int1d), particles_columns))
-- leaf (false), inner (false), idempotent (false)
where
  reads($particles#14023.particle_velocity), writes($particles#14023.particle_velocity), reads($particles#14023.position), writes($particles#14023.position), reads($particles#14023.position_ghost), reads($particles#14023.velocity_ghost), reads($particles#14023.velocity_t), writes($particles#14023.velocity_t), reads($particles#14023.velocity_t_ghost), reads($particles#14023.__valid), $dom#14021 <= $particles#14023
do
  for $p#14026 : int1d(particles_columns, $dom#14021) in $dom#14021 do
    if $particles#14023[$p#14026].__valid then
      $particles#14023[$p#14026].position = $particles#14023[$p#14026].position_ghost
      $particles#14023[$p#14026].particle_velocity = $particles#14023[$p#14026].velocity_ghost
      $particles#14023[$p#14026].velocity_t = $particles#14023[$p#14026].velocity_t_ghost
    else
    end
  end
end
task Particles_DeleteEscapingParticles($dom#14050 : region#3684(ispace#3684(int1d), particles_columns), $particles#14052 : region#3685(ispace#3685(int1d), particles_columns)) : int64
-- leaf (false), inner (false), idempotent (false)
where
  reads($particles#14052.position), reads($particles#14052.__valid), writes($particles#14052.__valid), reads($particles#14052.__valid), $dom#14050 <= $particles#14052
do
  var $acc#14092 : int64 = int64(0)
  for $p#14093 : int1d(particles_columns, $dom#14050) in $dom#14050 do
    if $particles#14052[$p#14093].__valid then
      var $min_x#14094 : double = double(0.6987052)
      var $max_x#14095 : double = (double(0.6987052)+double(0.5988125))
      var $min_y#14096 : double = double(0.99700342743682)
      var $max_y#14097 : double = (double(0.99700342743682)+double(0.84872084512635))
      var $min_z#14098 : double = double(0.70456500505464)
      var $max_z#14099 : double = (double(0.70456500505464)+double(0.83330648989071))
      var $pos#14100 : double[3] = $particles#14052[$p#14093].position
      if (((((($pos#14100[int32(0)]>$max_x#14095) or ($pos#14100[int32(0)]<$min_x#14094)) or ($pos#14100[int32(1)]>$max_y#14097)) or ($pos#14100[int32(1)]<$min_y#14096)) or ($pos#14100[int32(2)]>$max_z#14099)) or ($pos#14100[int32(2)]<$min_z#14098)) then
        $particles#14052[$p#14093].__valid = false
        $acc#14092 += (-int64(int32(1)))
      else
      end
    else
    end
  end
  return $acc#14092
end
task print______($14143 : double)
-- leaf (false), inner (false), idempotent (false)
  printf("\n Current time step: %2.6e s.\n", $14143)
end
task print_______($14149 : double, $14150 : double)
-- leaf (false), inner (false), idempotent (false)
  printf(" Min Flow Temp: %11.6f K. Max Flow Temp: %11.6f K.\n", $14149, $14150)
end
task print________($14157 : int64)
-- leaf (false), inner (false), idempotent (false)
  printf(" Current number of particles: %d.\n", $14157)
end
task print_________()
-- leaf (false), inner (false), idempotent (false)
  printf("\n")
end
task print__________()
-- leaf (false), inner (false), idempotent (false)
  printf("    Iter     Time(s)   Avg Press    Avg Temp      Avg KE  Particle T\n")
end
task print___________($14172 : int32, $14173 : double, $14174 : double, $14175 : double, $14176 : double, $14177 : double)
-- leaf (false), inner (false), idempotent (false)
  printf("%8d %11.6f %11.6f %11.6f %11.6f %11.6f\n", $14172, $14173, $14174, $14175, $14176, $14177)
end
terra Fluid_hdf5create_rho_pressure_velocity_(filename : &int8) : {}
    var fid : int32 = H5Fcreate(filename, [uint32](2), 0, 0)
    var dataSpace : int32
    var sizes : uint64[3]
    [&uint64](sizes)[2] = [uint64](378)
    [&uint64](sizes)[1] = [uint64](833)
    [&uint64](sizes)[0] = [uint64](734)
    dataSpace = H5Screate_simple(3, [&uint64](sizes), [&uint64](0))
    var hType : int32 = extern global H5T_IEEE_F64LE_g : int32
    var $dataSet : int32 = H5Dcreate2(fid, "rho", hType, dataSpace, 0, 0, 0)
    var hType$1 : int32 = extern global H5T_IEEE_F64LE_g : int32
    var $dataSet$1 : int32 = H5Dcreate2(fid, "pressure", hType$1, dataSpace, 0, 0, 0)
    var dims : uint64[1]
    [&uint64](dims)[0] = [uint64](3)
    var elemType : int32 = extern global H5T_IEEE_F64LE_g : int32
    var $arrayType : int32 = H5Tarray_create2(elemType, [uint32](1), [&uint64](dims))
    var hType$2 : int32 = $arrayType
    var $dataSet$2 : int32 = H5Dcreate2(fid, "velocity", hType$2, dataSpace, 0, 0, 0)
    H5Dclose($dataSet$2)
    H5Tclose($arrayType)
    H5Dclose($dataSet$1)
    H5Dclose($dataSet)
    H5Sclose(dataSpace)
    H5Fclose(fid)
end
terra particles_hdf5create_cell_position_particle_velocity_particle_temperature_diameter___valid_(filename : &int8) : {}
    var fid : int32 = H5Fcreate(filename, [uint32](2), 0, 0)
    var dataSpace : int32
    var sizes : uint64[1]
    [&uint64](sizes)[0] = [uint64](15 * 24)
    dataSpace = H5Screate_simple(1, [&uint64](sizes), [&uint64](0))
    var $int3dType : int32 = H5Tcreate([uint32](6), [uint64](24))
    H5Tinsert($int3dType, "x", [uint64](0), extern global H5T_STD_I64LE_g : int32)
    H5Tinsert($int3dType, "y", [uint64](8), extern global H5T_STD_I64LE_g : int32)
    H5Tinsert($int3dType, "z", [uint64](16), extern global H5T_STD_I64LE_g : int32)
    var $dataSet : int32 = H5Dcreate2(fid, "cell", $int3dType, dataSpace, 0, 0, 0)
    var dims : uint64[1]
    [&uint64](dims)[0] = [uint64](3)
    var elemType : int32 = extern global H5T_IEEE_F64LE_g : int32
    var $arrayType : int32 = H5Tarray_create2(elemType, [uint32](1), [&uint64](dims))
    var hType : int32 = $arrayType
    var $dataSet$1 : int32 = H5Dcreate2(fid, "position", hType, dataSpace, 0, 0, 0)
    var dims$1 : uint64[1]
    [&uint64](dims$1)[0] = [uint64](3)
    var elemType$1 : int32 = extern global H5T_IEEE_F64LE_g : int32
    var $arrayType$1 : int32 = H5Tarray_create2(elemType$1, [uint32](1), [&uint64](dims$1))
    var hType$1 : int32 = $arrayType$1
    var $dataSet$2 : int32 = H5Dcreate2(fid, "particle_velocity", hType$1, dataSpace, 0, 0, 0)
    var hType$2 : int32 = extern global H5T_IEEE_F64LE_g : int32
    var $dataSet$3 : int32 = H5Dcreate2(fid, "particle_temperature", hType$2, dataSpace, 0, 0, 0)
    var hType$3 : int32 = extern global H5T_IEEE_F64LE_g : int32
    var $dataSet$4 : int32 = H5Dcreate2(fid, "diameter", hType$3, dataSpace, 0, 0, 0)
    var hType$4 : int32 = extern global H5T_STD_U8LE_g : int32
    var $dataSet$5 : int32 = H5Dcreate2(fid, "__valid", hType$4, dataSpace, 0, 0, 0)
    H5Dclose($dataSet$5)
    H5Dclose($dataSet$4)
    H5Dclose($dataSet$3)
    H5Dclose($dataSet$2)
    H5Tclose($arrayType$1)
    H5Dclose($dataSet$1)
    H5Tclose($arrayType)
    H5Dclose($dataSet)
    H5Tclose($int3dType)
    H5Sclose(dataSpace)
    H5Fclose(fid)
end
task print____________($14286 : double)
-- leaf (false), inner (false), idempotent (false)
  printf("\n Current time step: %2.6e s.\n", $14286)
end
task print_____________($14292 : double, $14293 : double)
-- leaf (false), inner (false), idempotent (false)
  printf(" Min Flow Temp: %11.6f K. Max Flow Temp: %11.6f K.\n", $14292, $14293)
end
task print______________($14300 : int64)
-- leaf (false), inner (false), idempotent (false)
  printf(" Current number of particles: %d.\n", $14300)
end
task print_______________()
-- leaf (false), inner (false), idempotent (false)
  printf("\n")
end
task print________________()
-- leaf (false), inner (false), idempotent (false)
  printf("    Iter     Time(s)   Avg Press    Avg Temp      Avg KE  Particle T\n")
end
task print_________________($14315 : int32, $14316 : double, $14317 : double, $14318 : double, $14319 : double, $14320 : double)
-- leaf (false), inner (false), idempotent (false)
  printf("%8d %11.6f %11.6f %11.6f %11.6f %11.6f\n", $14315, $14316, $14317, $14318, $14319, $14320)
end
task main()
-- leaf (false), inner (false), idempotent (false)
  var $Flow_averageHeatSource#14334 : double = double(0)
  var $Flow_maxTemperature#14335 : double = double(0)
  var $Flow_averageDissipation#14336 : double = double(0)
  var $TimeIntegrator_timeOld#14337 : double = double(0)
  var $Flow_numberOfInteriorCells#14338 : int64 = int64(0)
  var $Particles_limit#14339 : int64 = int64(0)
  var $Flow_averageFe#14340 : double = double(0)
  var $Flow_areaInterior#14341 : double = double(0)
  var $Particles_averageTemperature#14342 : double = double(0)
  var $Flow_averagePD#14343 : double = double(0)
  var $maxH#14344 : double = double(0)
  var $maxV#14345 : double = double(0)
  var $TimeIntegrator_stage#14346 : int32 = int32(0)
  var $TimeIntegrator_timeStep#14347 : int32 = int32(0)
  var $maxC#14348 : double = double(0)
  var $Flow_averageK#14349 : double = double(0)
  var $Flow_minTemperature#14350 : double = double(0)
  var $Flow_averageKineticEnergy#14351 : double = double(0)
  var $Flow_averagePressure#14352 : double = double(0)
  var $TimeIntegrator_simTime#14353 : double = double(0)
  var $TimeIntegrator_deltaTime#14354 : double = double(0.0001)
  var $Flow_averageTemperature#14355 : double = double(0)
  var $Particles_number#14356 : int64 = int64(0)
  var $is#14357 : ispace#3729(int3d) = ispace(int3d, int3d({x = 378, y = 833, z = 734}))
  var $Fluid#14358 : region#3729(ispace#3729(int3d), Fluid_columns) = region($is#14357, Fluid_columns)
  var $Fluid_copy#14359 : region#3730(ispace#3729(int3d), Fluid_columns) = region($is#14357, Fluid_columns)
  var $colors#14360 : ispace#3730(int1d) = ispace(int1d, int1d(5))
  var $coloring#14361 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  legion_domain_point_coloring_color_domain($coloring#14361, legion_domain_point_t(int1d((1-1))), legion_domain_t(rect3d({lo = int3d({x = 0, y = 1, z = 733}), hi = int3d({x = 377, y = 831, z = 733})})))
  legion_domain_point_coloring_color_domain($coloring#14361, legion_domain_point_t(int1d((2-1))), legion_domain_t(rect3d({lo = int3d({x = 0, y = 0, z = 0}), hi = int3d({x = 377, y = 0, z = 733})})))
  legion_domain_point_coloring_color_domain($coloring#14361, legion_domain_point_t(int1d((3-1))), legion_domain_t(rect3d({lo = int3d({x = 0, y = 1, z = 1}), hi = int3d({x = 377, y = 831, z = 732})})))
  legion_domain_point_coloring_color_domain($coloring#14361, legion_domain_point_t(int1d((4-1))), legion_domain_t(rect3d({lo = int3d({x = 0, y = 1, z = 0}), hi = int3d({x = 377, y = 831, z = 0})})))
  legion_domain_point_coloring_color_domain($coloring#14361, legion_domain_point_t(int1d((5-1))), legion_domain_t(rect3d({lo = int3d({x = 0, y = 832, z = 0}), hi = int3d({x = 377, y = 832, z = 733})})))
  var $p#14362 : partition#3731(disjoint, $Fluid#14358, $colors#14360) = partition(disjoint, $Fluid#14358, $coloring#14361, $colors#14360)
  var $Fluid_boundary_zpos#14363 : region#3732(ispace#3731(int3d), Fluid_columns) = $p#14362[int1d((1-1))]
  var $Fluid_boundary_yneg#14364 : region#3733(ispace#3732(int3d), Fluid_columns) = $p#14362[int1d((2-1))]
  var $Fluid_interior#14365 : region#3734(ispace#3733(int3d), Fluid_columns) = $p#14362[int1d((3-1))]
  var $Fluid_boundary_zneg#14366 : region#3735(ispace#3734(int3d), Fluid_columns) = $p#14362[int1d((4-1))]
  var $Fluid_boundary_ypos#14367 : region#3736(ispace#3735(int3d), Fluid_columns) = $p#14362[int1d((5-1))]
  legion_domain_point_coloring_destroy($coloring#14361)
  var $is#14368 : ispace#3736(int1d) = ispace(int1d, int1d((15*24)))
  var $particles#14369 : region#3737(ispace#3736(int1d), particles_columns) = region($is#14368, particles_columns)
  var $particles_copy#14370 : region#3738(ispace#3736(int1d), particles_columns) = region($is#14368, particles_columns)
  var $primColors#14371 : ispace#3737(int3d) = ispace(int3d, int3d({2, 3, 4}))
  var $coloring#14372 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $c#14373 : int3d($primColors#14371) in $primColors#14371 do
    var $rect#14374 : rect3d = rect3d({lo = int3d({x = (0+(189*$c#14373.x)), y = (1+(277*$c#14373.y)), z = (1+(183*$c#14373.z))}), hi = int3d({x = ((0+(189*($c#14373.x+1)))-1), y = ((1+(277*($c#14373.y+1)))-1), z = ((1+(183*($c#14373.z+1)))-1)})})
    if ($c#14373.x==0) then
      $rect#14374.lo.x -= 0
    else
    end
    if ($c#14373.x==(2-1)) then
      $rect#14374.hi.x += 0
    else
    end
    if ($c#14373.y==0) then
      $rect#14374.lo.y -= 1
    else
    end
    if ($c#14373.y==(3-1)) then
      $rect#14374.hi.y += 1
    else
    end
    if ($c#14373.z==0) then
      $rect#14374.lo.z -= 1
    else
    end
    if ($c#14373.z==(4-1)) then
      $rect#14374.hi.z += 1
    else
    end
    legion_domain_point_coloring_color_domain($coloring#14372, legion_domain_point_t($c#14373), legion_domain_t($rect#14374))
  end
  var $Fluid_primPart#14375 : partition#3739(disjoint, $Fluid#14358, $primColors#14371) = partition(disjoint, $Fluid#14358, $coloring#14372, $primColors#14371)
  var $Fluid_copy_primPart#14376 : partition#3740(disjoint, $Fluid_copy#14359, $primColors#14371) = partition(disjoint, $Fluid_copy#14359, $coloring#14372, $primColors#14371)
  legion_domain_point_coloring_destroy($coloring#14372)
  var $coloring#14377 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14378 : int32 = 0, 4 do
    for $y#14379 : int32 = 0, 3 do
      for $x#14380 : int32 = 0, 2 do
        var $rBase#14381 : int64
        for $rStart#14382 : int1d(particles_columns, $particles#14369) in $particles#14369 do
          $rBase#14381 = int64(($rStart#14382+((((($z#14378*2)*3)+($y#14379*2))+$x#14380)*15)))
          break
        end
        legion_domain_point_coloring_color_domain($coloring#14377, legion_domain_point_t(int3d({$x#14380, $y#14379, $z#14378})), legion_domain_t(rect1d({$rBase#14381, (($rBase#14381+15)-1)})))
      end
    end
  end
  var $particles_primPart#14383 : partition#3741(disjoint, $particles#14369, $primColors#14371) = partition(disjoint, $particles#14369, $coloring#14377, $primColors#14371)
  var $particles_copy_primPart#14384 : partition#3742(disjoint, $particles_copy#14370, $primColors#14371) = partition(disjoint, $particles_copy#14370, $coloring#14377, $primColors#14371)
  legion_domain_point_coloring_destroy($coloring#14377)
  var $particles_queue_0#14385 : region#3743(ispace#3738(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14386 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14387 : int32 = 0, 4 do
    for $y#14388 : int32 = 0, 3 do
      for $x#14389 : int32 = 0, 2 do
        var $qBase#14390 : int64
        for $qStart#14391 : int1d(int8[376], $particles_queue_0#14385) in $particles_queue_0#14385 do
          $qBase#14390 = int64(($qStart#14391+((((($z#14387*2)*3)+($y#14388*2))+$x#14389)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14386, legion_domain_point_t(int3d({$x#14389, $y#14388, $z#14387})), legion_domain_t(rect1d({$qBase#14390, (($qBase#14390+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_0#14392 : partition#3744(disjoint, $particles_queue_0#14385, $primColors#14371) = partition(disjoint, $particles_queue_0#14385, $srcColoring#14386, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14386)
  var $dstColoring#14393 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14394 : int3d = int3d({0, 0, 1})
  for $c#14395 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14396 : int64
    for $qptr#14397 : int1d(int8[376], $15277) in $particles_qSrcPart_0#14392[((($c#14395-$colorOff#14394)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14396 = int64(int1d($qptr#14397))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14393, legion_domain_point_t($c#14395), legion_domain_t(rect1d({$srcBase#14396, (($srcBase#14396+1000)-1)})))
  end
  var $particles_qDstPart_0#14398 : partition#3746(aliased, $particles_queue_0#14385, $primColors#14371) = partition(aliased, $particles_queue_0#14385, $dstColoring#14393, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14393)
  var $particles_queue_1#14399 : region#3747(ispace#3740(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14400 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14401 : int32 = 0, 4 do
    for $y#14402 : int32 = 0, 3 do
      for $x#14403 : int32 = 0, 2 do
        var $qBase#14404 : int64
        for $qStart#14405 : int1d(int8[376], $particles_queue_1#14399) in $particles_queue_1#14399 do
          $qBase#14404 = int64(($qStart#14405+((((($z#14401*2)*3)+($y#14402*2))+$x#14403)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14400, legion_domain_point_t(int3d({$x#14403, $y#14402, $z#14401})), legion_domain_t(rect1d({$qBase#14404, (($qBase#14404+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_1#14406 : partition#3748(disjoint, $particles_queue_1#14399, $primColors#14371) = partition(disjoint, $particles_queue_1#14399, $srcColoring#14400, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14400)
  var $dstColoring#14407 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14408 : int3d = int3d({0, 0, -1})
  for $c#14409 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14410 : int64
    for $qptr#14411 : int1d(int8[376], $15291) in $particles_qSrcPart_1#14406[((($c#14409-$colorOff#14408)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14410 = int64(int1d($qptr#14411))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14407, legion_domain_point_t($c#14409), legion_domain_t(rect1d({$srcBase#14410, (($srcBase#14410+1000)-1)})))
  end
  var $particles_qDstPart_1#14412 : partition#3750(aliased, $particles_queue_1#14399, $primColors#14371) = partition(aliased, $particles_queue_1#14399, $dstColoring#14407, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14407)
  var $particles_queue_2#14413 : region#3751(ispace#3742(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14414 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14415 : int32 = 0, 4 do
    for $y#14416 : int32 = 0, 3 do
      for $x#14417 : int32 = 0, 2 do
        var $qBase#14418 : int64
        for $qStart#14419 : int1d(int8[376], $particles_queue_2#14413) in $particles_queue_2#14413 do
          $qBase#14418 = int64(($qStart#14419+((((($z#14415*2)*3)+($y#14416*2))+$x#14417)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14414, legion_domain_point_t(int3d({$x#14417, $y#14416, $z#14415})), legion_domain_t(rect1d({$qBase#14418, (($qBase#14418+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_2#14420 : partition#3752(disjoint, $particles_queue_2#14413, $primColors#14371) = partition(disjoint, $particles_queue_2#14413, $srcColoring#14414, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14414)
  var $dstColoring#14421 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14422 : int3d = int3d({0, 1, 0})
  for $c#14423 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14424 : int64
    for $qptr#14425 : int1d(int8[376], $15305) in $particles_qSrcPart_2#14420[((($c#14423-$colorOff#14422)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14424 = int64(int1d($qptr#14425))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14421, legion_domain_point_t($c#14423), legion_domain_t(rect1d({$srcBase#14424, (($srcBase#14424+1000)-1)})))
  end
  var $particles_qDstPart_2#14426 : partition#3754(aliased, $particles_queue_2#14413, $primColors#14371) = partition(aliased, $particles_queue_2#14413, $dstColoring#14421, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14421)
  var $particles_queue_3#14427 : region#3755(ispace#3744(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14428 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14429 : int32 = 0, 4 do
    for $y#14430 : int32 = 0, 3 do
      for $x#14431 : int32 = 0, 2 do
        var $qBase#14432 : int64
        for $qStart#14433 : int1d(int8[376], $particles_queue_3#14427) in $particles_queue_3#14427 do
          $qBase#14432 = int64(($qStart#14433+((((($z#14429*2)*3)+($y#14430*2))+$x#14431)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14428, legion_domain_point_t(int3d({$x#14431, $y#14430, $z#14429})), legion_domain_t(rect1d({$qBase#14432, (($qBase#14432+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_3#14434 : partition#3756(disjoint, $particles_queue_3#14427, $primColors#14371) = partition(disjoint, $particles_queue_3#14427, $srcColoring#14428, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14428)
  var $dstColoring#14435 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14436 : int3d = int3d({0, 1, 1})
  for $c#14437 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14438 : int64
    for $qptr#14439 : int1d(int8[376], $15319) in $particles_qSrcPart_3#14434[((($c#14437-$colorOff#14436)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14438 = int64(int1d($qptr#14439))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14435, legion_domain_point_t($c#14437), legion_domain_t(rect1d({$srcBase#14438, (($srcBase#14438+1000)-1)})))
  end
  var $particles_qDstPart_3#14440 : partition#3758(aliased, $particles_queue_3#14427, $primColors#14371) = partition(aliased, $particles_queue_3#14427, $dstColoring#14435, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14435)
  var $particles_queue_4#14441 : region#3759(ispace#3746(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14442 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14443 : int32 = 0, 4 do
    for $y#14444 : int32 = 0, 3 do
      for $x#14445 : int32 = 0, 2 do
        var $qBase#14446 : int64
        for $qStart#14447 : int1d(int8[376], $particles_queue_4#14441) in $particles_queue_4#14441 do
          $qBase#14446 = int64(($qStart#14447+((((($z#14443*2)*3)+($y#14444*2))+$x#14445)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14442, legion_domain_point_t(int3d({$x#14445, $y#14444, $z#14443})), legion_domain_t(rect1d({$qBase#14446, (($qBase#14446+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_4#14448 : partition#3760(disjoint, $particles_queue_4#14441, $primColors#14371) = partition(disjoint, $particles_queue_4#14441, $srcColoring#14442, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14442)
  var $dstColoring#14449 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14450 : int3d = int3d({0, 1, -1})
  for $c#14451 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14452 : int64
    for $qptr#14453 : int1d(int8[376], $15333) in $particles_qSrcPart_4#14448[((($c#14451-$colorOff#14450)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14452 = int64(int1d($qptr#14453))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14449, legion_domain_point_t($c#14451), legion_domain_t(rect1d({$srcBase#14452, (($srcBase#14452+1000)-1)})))
  end
  var $particles_qDstPart_4#14454 : partition#3762(aliased, $particles_queue_4#14441, $primColors#14371) = partition(aliased, $particles_queue_4#14441, $dstColoring#14449, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14449)
  var $particles_queue_5#14455 : region#3763(ispace#3748(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14456 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14457 : int32 = 0, 4 do
    for $y#14458 : int32 = 0, 3 do
      for $x#14459 : int32 = 0, 2 do
        var $qBase#14460 : int64
        for $qStart#14461 : int1d(int8[376], $particles_queue_5#14455) in $particles_queue_5#14455 do
          $qBase#14460 = int64(($qStart#14461+((((($z#14457*2)*3)+($y#14458*2))+$x#14459)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14456, legion_domain_point_t(int3d({$x#14459, $y#14458, $z#14457})), legion_domain_t(rect1d({$qBase#14460, (($qBase#14460+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_5#14462 : partition#3764(disjoint, $particles_queue_5#14455, $primColors#14371) = partition(disjoint, $particles_queue_5#14455, $srcColoring#14456, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14456)
  var $dstColoring#14463 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14464 : int3d = int3d({0, -1, 0})
  for $c#14465 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14466 : int64
    for $qptr#14467 : int1d(int8[376], $15347) in $particles_qSrcPart_5#14462[((($c#14465-$colorOff#14464)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14466 = int64(int1d($qptr#14467))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14463, legion_domain_point_t($c#14465), legion_domain_t(rect1d({$srcBase#14466, (($srcBase#14466+1000)-1)})))
  end
  var $particles_qDstPart_5#14468 : partition#3766(aliased, $particles_queue_5#14455, $primColors#14371) = partition(aliased, $particles_queue_5#14455, $dstColoring#14463, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14463)
  var $particles_queue_6#14469 : region#3767(ispace#3750(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14470 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14471 : int32 = 0, 4 do
    for $y#14472 : int32 = 0, 3 do
      for $x#14473 : int32 = 0, 2 do
        var $qBase#14474 : int64
        for $qStart#14475 : int1d(int8[376], $particles_queue_6#14469) in $particles_queue_6#14469 do
          $qBase#14474 = int64(($qStart#14475+((((($z#14471*2)*3)+($y#14472*2))+$x#14473)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14470, legion_domain_point_t(int3d({$x#14473, $y#14472, $z#14471})), legion_domain_t(rect1d({$qBase#14474, (($qBase#14474+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_6#14476 : partition#3768(disjoint, $particles_queue_6#14469, $primColors#14371) = partition(disjoint, $particles_queue_6#14469, $srcColoring#14470, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14470)
  var $dstColoring#14477 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14478 : int3d = int3d({0, -1, 1})
  for $c#14479 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14480 : int64
    for $qptr#14481 : int1d(int8[376], $15361) in $particles_qSrcPart_6#14476[((($c#14479-$colorOff#14478)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14480 = int64(int1d($qptr#14481))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14477, legion_domain_point_t($c#14479), legion_domain_t(rect1d({$srcBase#14480, (($srcBase#14480+1000)-1)})))
  end
  var $particles_qDstPart_6#14482 : partition#3770(aliased, $particles_queue_6#14469, $primColors#14371) = partition(aliased, $particles_queue_6#14469, $dstColoring#14477, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14477)
  var $particles_queue_7#14483 : region#3771(ispace#3752(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14484 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14485 : int32 = 0, 4 do
    for $y#14486 : int32 = 0, 3 do
      for $x#14487 : int32 = 0, 2 do
        var $qBase#14488 : int64
        for $qStart#14489 : int1d(int8[376], $particles_queue_7#14483) in $particles_queue_7#14483 do
          $qBase#14488 = int64(($qStart#14489+((((($z#14485*2)*3)+($y#14486*2))+$x#14487)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14484, legion_domain_point_t(int3d({$x#14487, $y#14486, $z#14485})), legion_domain_t(rect1d({$qBase#14488, (($qBase#14488+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_7#14490 : partition#3772(disjoint, $particles_queue_7#14483, $primColors#14371) = partition(disjoint, $particles_queue_7#14483, $srcColoring#14484, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14484)
  var $dstColoring#14491 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14492 : int3d = int3d({0, -1, -1})
  for $c#14493 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14494 : int64
    for $qptr#14495 : int1d(int8[376], $15375) in $particles_qSrcPart_7#14490[((($c#14493-$colorOff#14492)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14494 = int64(int1d($qptr#14495))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14491, legion_domain_point_t($c#14493), legion_domain_t(rect1d({$srcBase#14494, (($srcBase#14494+1000)-1)})))
  end
  var $particles_qDstPart_7#14496 : partition#3774(aliased, $particles_queue_7#14483, $primColors#14371) = partition(aliased, $particles_queue_7#14483, $dstColoring#14491, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14491)
  var $particles_queue_8#14497 : region#3775(ispace#3754(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14498 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14499 : int32 = 0, 4 do
    for $y#14500 : int32 = 0, 3 do
      for $x#14501 : int32 = 0, 2 do
        var $qBase#14502 : int64
        for $qStart#14503 : int1d(int8[376], $particles_queue_8#14497) in $particles_queue_8#14497 do
          $qBase#14502 = int64(($qStart#14503+((((($z#14499*2)*3)+($y#14500*2))+$x#14501)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14498, legion_domain_point_t(int3d({$x#14501, $y#14500, $z#14499})), legion_domain_t(rect1d({$qBase#14502, (($qBase#14502+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_8#14504 : partition#3776(disjoint, $particles_queue_8#14497, $primColors#14371) = partition(disjoint, $particles_queue_8#14497, $srcColoring#14498, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14498)
  var $dstColoring#14505 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14506 : int3d = int3d({1, 0, 0})
  for $c#14507 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14508 : int64
    for $qptr#14509 : int1d(int8[376], $15389) in $particles_qSrcPart_8#14504[((($c#14507-$colorOff#14506)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14508 = int64(int1d($qptr#14509))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14505, legion_domain_point_t($c#14507), legion_domain_t(rect1d({$srcBase#14508, (($srcBase#14508+1000)-1)})))
  end
  var $particles_qDstPart_8#14510 : partition#3778(aliased, $particles_queue_8#14497, $primColors#14371) = partition(aliased, $particles_queue_8#14497, $dstColoring#14505, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14505)
  var $particles_queue_9#14511 : region#3779(ispace#3756(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14512 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14513 : int32 = 0, 4 do
    for $y#14514 : int32 = 0, 3 do
      for $x#14515 : int32 = 0, 2 do
        var $qBase#14516 : int64
        for $qStart#14517 : int1d(int8[376], $particles_queue_9#14511) in $particles_queue_9#14511 do
          $qBase#14516 = int64(($qStart#14517+((((($z#14513*2)*3)+($y#14514*2))+$x#14515)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14512, legion_domain_point_t(int3d({$x#14515, $y#14514, $z#14513})), legion_domain_t(rect1d({$qBase#14516, (($qBase#14516+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_9#14518 : partition#3780(disjoint, $particles_queue_9#14511, $primColors#14371) = partition(disjoint, $particles_queue_9#14511, $srcColoring#14512, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14512)
  var $dstColoring#14519 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14520 : int3d = int3d({1, 0, 1})
  for $c#14521 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14522 : int64
    for $qptr#14523 : int1d(int8[376], $15403) in $particles_qSrcPart_9#14518[((($c#14521-$colorOff#14520)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14522 = int64(int1d($qptr#14523))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14519, legion_domain_point_t($c#14521), legion_domain_t(rect1d({$srcBase#14522, (($srcBase#14522+1000)-1)})))
  end
  var $particles_qDstPart_9#14524 : partition#3782(aliased, $particles_queue_9#14511, $primColors#14371) = partition(aliased, $particles_queue_9#14511, $dstColoring#14519, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14519)
  var $particles_queue_10#14525 : region#3783(ispace#3758(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14526 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14527 : int32 = 0, 4 do
    for $y#14528 : int32 = 0, 3 do
      for $x#14529 : int32 = 0, 2 do
        var $qBase#14530 : int64
        for $qStart#14531 : int1d(int8[376], $particles_queue_10#14525) in $particles_queue_10#14525 do
          $qBase#14530 = int64(($qStart#14531+((((($z#14527*2)*3)+($y#14528*2))+$x#14529)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14526, legion_domain_point_t(int3d({$x#14529, $y#14528, $z#14527})), legion_domain_t(rect1d({$qBase#14530, (($qBase#14530+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_10#14532 : partition#3784(disjoint, $particles_queue_10#14525, $primColors#14371) = partition(disjoint, $particles_queue_10#14525, $srcColoring#14526, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14526)
  var $dstColoring#14533 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14534 : int3d = int3d({1, 0, -1})
  for $c#14535 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14536 : int64
    for $qptr#14537 : int1d(int8[376], $15417) in $particles_qSrcPart_10#14532[((($c#14535-$colorOff#14534)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14536 = int64(int1d($qptr#14537))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14533, legion_domain_point_t($c#14535), legion_domain_t(rect1d({$srcBase#14536, (($srcBase#14536+1000)-1)})))
  end
  var $particles_qDstPart_10#14538 : partition#3786(aliased, $particles_queue_10#14525, $primColors#14371) = partition(aliased, $particles_queue_10#14525, $dstColoring#14533, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14533)
  var $particles_queue_11#14539 : region#3787(ispace#3760(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14540 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14541 : int32 = 0, 4 do
    for $y#14542 : int32 = 0, 3 do
      for $x#14543 : int32 = 0, 2 do
        var $qBase#14544 : int64
        for $qStart#14545 : int1d(int8[376], $particles_queue_11#14539) in $particles_queue_11#14539 do
          $qBase#14544 = int64(($qStart#14545+((((($z#14541*2)*3)+($y#14542*2))+$x#14543)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14540, legion_domain_point_t(int3d({$x#14543, $y#14542, $z#14541})), legion_domain_t(rect1d({$qBase#14544, (($qBase#14544+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_11#14546 : partition#3788(disjoint, $particles_queue_11#14539, $primColors#14371) = partition(disjoint, $particles_queue_11#14539, $srcColoring#14540, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14540)
  var $dstColoring#14547 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14548 : int3d = int3d({1, 1, 0})
  for $c#14549 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14550 : int64
    for $qptr#14551 : int1d(int8[376], $15431) in $particles_qSrcPart_11#14546[((($c#14549-$colorOff#14548)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14550 = int64(int1d($qptr#14551))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14547, legion_domain_point_t($c#14549), legion_domain_t(rect1d({$srcBase#14550, (($srcBase#14550+1000)-1)})))
  end
  var $particles_qDstPart_11#14552 : partition#3790(aliased, $particles_queue_11#14539, $primColors#14371) = partition(aliased, $particles_queue_11#14539, $dstColoring#14547, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14547)
  var $particles_queue_12#14553 : region#3791(ispace#3762(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14554 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14555 : int32 = 0, 4 do
    for $y#14556 : int32 = 0, 3 do
      for $x#14557 : int32 = 0, 2 do
        var $qBase#14558 : int64
        for $qStart#14559 : int1d(int8[376], $particles_queue_12#14553) in $particles_queue_12#14553 do
          $qBase#14558 = int64(($qStart#14559+((((($z#14555*2)*3)+($y#14556*2))+$x#14557)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14554, legion_domain_point_t(int3d({$x#14557, $y#14556, $z#14555})), legion_domain_t(rect1d({$qBase#14558, (($qBase#14558+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_12#14560 : partition#3792(disjoint, $particles_queue_12#14553, $primColors#14371) = partition(disjoint, $particles_queue_12#14553, $srcColoring#14554, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14554)
  var $dstColoring#14561 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14562 : int3d = int3d({1, 1, 1})
  for $c#14563 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14564 : int64
    for $qptr#14565 : int1d(int8[376], $15445) in $particles_qSrcPart_12#14560[((($c#14563-$colorOff#14562)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14564 = int64(int1d($qptr#14565))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14561, legion_domain_point_t($c#14563), legion_domain_t(rect1d({$srcBase#14564, (($srcBase#14564+1000)-1)})))
  end
  var $particles_qDstPart_12#14566 : partition#3794(aliased, $particles_queue_12#14553, $primColors#14371) = partition(aliased, $particles_queue_12#14553, $dstColoring#14561, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14561)
  var $particles_queue_13#14567 : region#3795(ispace#3764(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14568 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14569 : int32 = 0, 4 do
    for $y#14570 : int32 = 0, 3 do
      for $x#14571 : int32 = 0, 2 do
        var $qBase#14572 : int64
        for $qStart#14573 : int1d(int8[376], $particles_queue_13#14567) in $particles_queue_13#14567 do
          $qBase#14572 = int64(($qStart#14573+((((($z#14569*2)*3)+($y#14570*2))+$x#14571)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14568, legion_domain_point_t(int3d({$x#14571, $y#14570, $z#14569})), legion_domain_t(rect1d({$qBase#14572, (($qBase#14572+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_13#14574 : partition#3796(disjoint, $particles_queue_13#14567, $primColors#14371) = partition(disjoint, $particles_queue_13#14567, $srcColoring#14568, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14568)
  var $dstColoring#14575 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14576 : int3d = int3d({1, 1, -1})
  for $c#14577 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14578 : int64
    for $qptr#14579 : int1d(int8[376], $15459) in $particles_qSrcPart_13#14574[((($c#14577-$colorOff#14576)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14578 = int64(int1d($qptr#14579))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14575, legion_domain_point_t($c#14577), legion_domain_t(rect1d({$srcBase#14578, (($srcBase#14578+1000)-1)})))
  end
  var $particles_qDstPart_13#14580 : partition#3798(aliased, $particles_queue_13#14567, $primColors#14371) = partition(aliased, $particles_queue_13#14567, $dstColoring#14575, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14575)
  var $particles_queue_14#14581 : region#3799(ispace#3766(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14582 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14583 : int32 = 0, 4 do
    for $y#14584 : int32 = 0, 3 do
      for $x#14585 : int32 = 0, 2 do
        var $qBase#14586 : int64
        for $qStart#14587 : int1d(int8[376], $particles_queue_14#14581) in $particles_queue_14#14581 do
          $qBase#14586 = int64(($qStart#14587+((((($z#14583*2)*3)+($y#14584*2))+$x#14585)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14582, legion_domain_point_t(int3d({$x#14585, $y#14584, $z#14583})), legion_domain_t(rect1d({$qBase#14586, (($qBase#14586+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_14#14588 : partition#3800(disjoint, $particles_queue_14#14581, $primColors#14371) = partition(disjoint, $particles_queue_14#14581, $srcColoring#14582, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14582)
  var $dstColoring#14589 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14590 : int3d = int3d({1, -1, 0})
  for $c#14591 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14592 : int64
    for $qptr#14593 : int1d(int8[376], $15473) in $particles_qSrcPart_14#14588[((($c#14591-$colorOff#14590)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14592 = int64(int1d($qptr#14593))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14589, legion_domain_point_t($c#14591), legion_domain_t(rect1d({$srcBase#14592, (($srcBase#14592+1000)-1)})))
  end
  var $particles_qDstPart_14#14594 : partition#3802(aliased, $particles_queue_14#14581, $primColors#14371) = partition(aliased, $particles_queue_14#14581, $dstColoring#14589, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14589)
  var $particles_queue_15#14595 : region#3803(ispace#3768(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14596 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14597 : int32 = 0, 4 do
    for $y#14598 : int32 = 0, 3 do
      for $x#14599 : int32 = 0, 2 do
        var $qBase#14600 : int64
        for $qStart#14601 : int1d(int8[376], $particles_queue_15#14595) in $particles_queue_15#14595 do
          $qBase#14600 = int64(($qStart#14601+((((($z#14597*2)*3)+($y#14598*2))+$x#14599)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14596, legion_domain_point_t(int3d({$x#14599, $y#14598, $z#14597})), legion_domain_t(rect1d({$qBase#14600, (($qBase#14600+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_15#14602 : partition#3804(disjoint, $particles_queue_15#14595, $primColors#14371) = partition(disjoint, $particles_queue_15#14595, $srcColoring#14596, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14596)
  var $dstColoring#14603 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14604 : int3d = int3d({1, -1, 1})
  for $c#14605 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14606 : int64
    for $qptr#14607 : int1d(int8[376], $15487) in $particles_qSrcPart_15#14602[((($c#14605-$colorOff#14604)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14606 = int64(int1d($qptr#14607))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14603, legion_domain_point_t($c#14605), legion_domain_t(rect1d({$srcBase#14606, (($srcBase#14606+1000)-1)})))
  end
  var $particles_qDstPart_15#14608 : partition#3806(aliased, $particles_queue_15#14595, $primColors#14371) = partition(aliased, $particles_queue_15#14595, $dstColoring#14603, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14603)
  var $particles_queue_16#14609 : region#3807(ispace#3770(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14610 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14611 : int32 = 0, 4 do
    for $y#14612 : int32 = 0, 3 do
      for $x#14613 : int32 = 0, 2 do
        var $qBase#14614 : int64
        for $qStart#14615 : int1d(int8[376], $particles_queue_16#14609) in $particles_queue_16#14609 do
          $qBase#14614 = int64(($qStart#14615+((((($z#14611*2)*3)+($y#14612*2))+$x#14613)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14610, legion_domain_point_t(int3d({$x#14613, $y#14612, $z#14611})), legion_domain_t(rect1d({$qBase#14614, (($qBase#14614+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_16#14616 : partition#3808(disjoint, $particles_queue_16#14609, $primColors#14371) = partition(disjoint, $particles_queue_16#14609, $srcColoring#14610, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14610)
  var $dstColoring#14617 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14618 : int3d = int3d({1, -1, -1})
  for $c#14619 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14620 : int64
    for $qptr#14621 : int1d(int8[376], $15501) in $particles_qSrcPart_16#14616[((($c#14619-$colorOff#14618)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14620 = int64(int1d($qptr#14621))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14617, legion_domain_point_t($c#14619), legion_domain_t(rect1d({$srcBase#14620, (($srcBase#14620+1000)-1)})))
  end
  var $particles_qDstPart_16#14622 : partition#3810(aliased, $particles_queue_16#14609, $primColors#14371) = partition(aliased, $particles_queue_16#14609, $dstColoring#14617, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14617)
  var $particles_queue_17#14623 : region#3811(ispace#3772(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14624 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14625 : int32 = 0, 4 do
    for $y#14626 : int32 = 0, 3 do
      for $x#14627 : int32 = 0, 2 do
        var $qBase#14628 : int64
        for $qStart#14629 : int1d(int8[376], $particles_queue_17#14623) in $particles_queue_17#14623 do
          $qBase#14628 = int64(($qStart#14629+((((($z#14625*2)*3)+($y#14626*2))+$x#14627)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14624, legion_domain_point_t(int3d({$x#14627, $y#14626, $z#14625})), legion_domain_t(rect1d({$qBase#14628, (($qBase#14628+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_17#14630 : partition#3812(disjoint, $particles_queue_17#14623, $primColors#14371) = partition(disjoint, $particles_queue_17#14623, $srcColoring#14624, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14624)
  var $dstColoring#14631 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14632 : int3d = int3d({-1, 0, 0})
  for $c#14633 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14634 : int64
    for $qptr#14635 : int1d(int8[376], $15515) in $particles_qSrcPart_17#14630[((($c#14633-$colorOff#14632)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14634 = int64(int1d($qptr#14635))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14631, legion_domain_point_t($c#14633), legion_domain_t(rect1d({$srcBase#14634, (($srcBase#14634+1000)-1)})))
  end
  var $particles_qDstPart_17#14636 : partition#3814(aliased, $particles_queue_17#14623, $primColors#14371) = partition(aliased, $particles_queue_17#14623, $dstColoring#14631, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14631)
  var $particles_queue_18#14637 : region#3815(ispace#3774(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14638 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14639 : int32 = 0, 4 do
    for $y#14640 : int32 = 0, 3 do
      for $x#14641 : int32 = 0, 2 do
        var $qBase#14642 : int64
        for $qStart#14643 : int1d(int8[376], $particles_queue_18#14637) in $particles_queue_18#14637 do
          $qBase#14642 = int64(($qStart#14643+((((($z#14639*2)*3)+($y#14640*2))+$x#14641)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14638, legion_domain_point_t(int3d({$x#14641, $y#14640, $z#14639})), legion_domain_t(rect1d({$qBase#14642, (($qBase#14642+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_18#14644 : partition#3816(disjoint, $particles_queue_18#14637, $primColors#14371) = partition(disjoint, $particles_queue_18#14637, $srcColoring#14638, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14638)
  var $dstColoring#14645 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14646 : int3d = int3d({-1, 0, 1})
  for $c#14647 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14648 : int64
    for $qptr#14649 : int1d(int8[376], $15529) in $particles_qSrcPart_18#14644[((($c#14647-$colorOff#14646)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14648 = int64(int1d($qptr#14649))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14645, legion_domain_point_t($c#14647), legion_domain_t(rect1d({$srcBase#14648, (($srcBase#14648+1000)-1)})))
  end
  var $particles_qDstPart_18#14650 : partition#3818(aliased, $particles_queue_18#14637, $primColors#14371) = partition(aliased, $particles_queue_18#14637, $dstColoring#14645, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14645)
  var $particles_queue_19#14651 : region#3819(ispace#3776(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14652 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14653 : int32 = 0, 4 do
    for $y#14654 : int32 = 0, 3 do
      for $x#14655 : int32 = 0, 2 do
        var $qBase#14656 : int64
        for $qStart#14657 : int1d(int8[376], $particles_queue_19#14651) in $particles_queue_19#14651 do
          $qBase#14656 = int64(($qStart#14657+((((($z#14653*2)*3)+($y#14654*2))+$x#14655)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14652, legion_domain_point_t(int3d({$x#14655, $y#14654, $z#14653})), legion_domain_t(rect1d({$qBase#14656, (($qBase#14656+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_19#14658 : partition#3820(disjoint, $particles_queue_19#14651, $primColors#14371) = partition(disjoint, $particles_queue_19#14651, $srcColoring#14652, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14652)
  var $dstColoring#14659 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14660 : int3d = int3d({-1, 0, -1})
  for $c#14661 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14662 : int64
    for $qptr#14663 : int1d(int8[376], $15543) in $particles_qSrcPart_19#14658[((($c#14661-$colorOff#14660)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14662 = int64(int1d($qptr#14663))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14659, legion_domain_point_t($c#14661), legion_domain_t(rect1d({$srcBase#14662, (($srcBase#14662+1000)-1)})))
  end
  var $particles_qDstPart_19#14664 : partition#3822(aliased, $particles_queue_19#14651, $primColors#14371) = partition(aliased, $particles_queue_19#14651, $dstColoring#14659, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14659)
  var $particles_queue_20#14665 : region#3823(ispace#3778(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14666 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14667 : int32 = 0, 4 do
    for $y#14668 : int32 = 0, 3 do
      for $x#14669 : int32 = 0, 2 do
        var $qBase#14670 : int64
        for $qStart#14671 : int1d(int8[376], $particles_queue_20#14665) in $particles_queue_20#14665 do
          $qBase#14670 = int64(($qStart#14671+((((($z#14667*2)*3)+($y#14668*2))+$x#14669)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14666, legion_domain_point_t(int3d({$x#14669, $y#14668, $z#14667})), legion_domain_t(rect1d({$qBase#14670, (($qBase#14670+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_20#14672 : partition#3824(disjoint, $particles_queue_20#14665, $primColors#14371) = partition(disjoint, $particles_queue_20#14665, $srcColoring#14666, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14666)
  var $dstColoring#14673 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14674 : int3d = int3d({-1, 1, 0})
  for $c#14675 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14676 : int64
    for $qptr#14677 : int1d(int8[376], $15557) in $particles_qSrcPart_20#14672[((($c#14675-$colorOff#14674)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14676 = int64(int1d($qptr#14677))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14673, legion_domain_point_t($c#14675), legion_domain_t(rect1d({$srcBase#14676, (($srcBase#14676+1000)-1)})))
  end
  var $particles_qDstPart_20#14678 : partition#3826(aliased, $particles_queue_20#14665, $primColors#14371) = partition(aliased, $particles_queue_20#14665, $dstColoring#14673, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14673)
  var $particles_queue_21#14679 : region#3827(ispace#3780(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14680 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14681 : int32 = 0, 4 do
    for $y#14682 : int32 = 0, 3 do
      for $x#14683 : int32 = 0, 2 do
        var $qBase#14684 : int64
        for $qStart#14685 : int1d(int8[376], $particles_queue_21#14679) in $particles_queue_21#14679 do
          $qBase#14684 = int64(($qStart#14685+((((($z#14681*2)*3)+($y#14682*2))+$x#14683)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14680, legion_domain_point_t(int3d({$x#14683, $y#14682, $z#14681})), legion_domain_t(rect1d({$qBase#14684, (($qBase#14684+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_21#14686 : partition#3828(disjoint, $particles_queue_21#14679, $primColors#14371) = partition(disjoint, $particles_queue_21#14679, $srcColoring#14680, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14680)
  var $dstColoring#14687 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14688 : int3d = int3d({-1, 1, 1})
  for $c#14689 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14690 : int64
    for $qptr#14691 : int1d(int8[376], $15571) in $particles_qSrcPart_21#14686[((($c#14689-$colorOff#14688)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14690 = int64(int1d($qptr#14691))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14687, legion_domain_point_t($c#14689), legion_domain_t(rect1d({$srcBase#14690, (($srcBase#14690+1000)-1)})))
  end
  var $particles_qDstPart_21#14692 : partition#3830(aliased, $particles_queue_21#14679, $primColors#14371) = partition(aliased, $particles_queue_21#14679, $dstColoring#14687, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14687)
  var $particles_queue_22#14693 : region#3831(ispace#3782(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14694 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14695 : int32 = 0, 4 do
    for $y#14696 : int32 = 0, 3 do
      for $x#14697 : int32 = 0, 2 do
        var $qBase#14698 : int64
        for $qStart#14699 : int1d(int8[376], $particles_queue_22#14693) in $particles_queue_22#14693 do
          $qBase#14698 = int64(($qStart#14699+((((($z#14695*2)*3)+($y#14696*2))+$x#14697)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14694, legion_domain_point_t(int3d({$x#14697, $y#14696, $z#14695})), legion_domain_t(rect1d({$qBase#14698, (($qBase#14698+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_22#14700 : partition#3832(disjoint, $particles_queue_22#14693, $primColors#14371) = partition(disjoint, $particles_queue_22#14693, $srcColoring#14694, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14694)
  var $dstColoring#14701 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14702 : int3d = int3d({-1, 1, -1})
  for $c#14703 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14704 : int64
    for $qptr#14705 : int1d(int8[376], $15585) in $particles_qSrcPart_22#14700[((($c#14703-$colorOff#14702)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14704 = int64(int1d($qptr#14705))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14701, legion_domain_point_t($c#14703), legion_domain_t(rect1d({$srcBase#14704, (($srcBase#14704+1000)-1)})))
  end
  var $particles_qDstPart_22#14706 : partition#3834(aliased, $particles_queue_22#14693, $primColors#14371) = partition(aliased, $particles_queue_22#14693, $dstColoring#14701, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14701)
  var $particles_queue_23#14707 : region#3835(ispace#3784(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14708 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14709 : int32 = 0, 4 do
    for $y#14710 : int32 = 0, 3 do
      for $x#14711 : int32 = 0, 2 do
        var $qBase#14712 : int64
        for $qStart#14713 : int1d(int8[376], $particles_queue_23#14707) in $particles_queue_23#14707 do
          $qBase#14712 = int64(($qStart#14713+((((($z#14709*2)*3)+($y#14710*2))+$x#14711)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14708, legion_domain_point_t(int3d({$x#14711, $y#14710, $z#14709})), legion_domain_t(rect1d({$qBase#14712, (($qBase#14712+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_23#14714 : partition#3836(disjoint, $particles_queue_23#14707, $primColors#14371) = partition(disjoint, $particles_queue_23#14707, $srcColoring#14708, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14708)
  var $dstColoring#14715 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14716 : int3d = int3d({-1, -1, 0})
  for $c#14717 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14718 : int64
    for $qptr#14719 : int1d(int8[376], $15599) in $particles_qSrcPart_23#14714[((($c#14717-$colorOff#14716)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14718 = int64(int1d($qptr#14719))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14715, legion_domain_point_t($c#14717), legion_domain_t(rect1d({$srcBase#14718, (($srcBase#14718+1000)-1)})))
  end
  var $particles_qDstPart_23#14720 : partition#3838(aliased, $particles_queue_23#14707, $primColors#14371) = partition(aliased, $particles_queue_23#14707, $dstColoring#14715, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14715)
  var $particles_queue_24#14721 : region#3839(ispace#3786(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14722 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14723 : int32 = 0, 4 do
    for $y#14724 : int32 = 0, 3 do
      for $x#14725 : int32 = 0, 2 do
        var $qBase#14726 : int64
        for $qStart#14727 : int1d(int8[376], $particles_queue_24#14721) in $particles_queue_24#14721 do
          $qBase#14726 = int64(($qStart#14727+((((($z#14723*2)*3)+($y#14724*2))+$x#14725)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14722, legion_domain_point_t(int3d({$x#14725, $y#14724, $z#14723})), legion_domain_t(rect1d({$qBase#14726, (($qBase#14726+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_24#14728 : partition#3840(disjoint, $particles_queue_24#14721, $primColors#14371) = partition(disjoint, $particles_queue_24#14721, $srcColoring#14722, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14722)
  var $dstColoring#14729 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14730 : int3d = int3d({-1, -1, 1})
  for $c#14731 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14732 : int64
    for $qptr#14733 : int1d(int8[376], $15613) in $particles_qSrcPart_24#14728[((($c#14731-$colorOff#14730)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14732 = int64(int1d($qptr#14733))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14729, legion_domain_point_t($c#14731), legion_domain_t(rect1d({$srcBase#14732, (($srcBase#14732+1000)-1)})))
  end
  var $particles_qDstPart_24#14734 : partition#3842(aliased, $particles_queue_24#14721, $primColors#14371) = partition(aliased, $particles_queue_24#14721, $dstColoring#14729, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14729)
  var $particles_queue_25#14735 : region#3843(ispace#3788(int1d), int8[376]) = region(ispace(int1d, int1d(24000)), int8[376])
  var $srcColoring#14736 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  for $z#14737 : int32 = 0, 4 do
    for $y#14738 : int32 = 0, 3 do
      for $x#14739 : int32 = 0, 2 do
        var $qBase#14740 : int64
        for $qStart#14741 : int1d(int8[376], $particles_queue_25#14735) in $particles_queue_25#14735 do
          $qBase#14740 = int64(($qStart#14741+((((($z#14737*2)*3)+($y#14738*2))+$x#14739)*1000)))
          break
        end
        legion_domain_point_coloring_color_domain($srcColoring#14736, legion_domain_point_t(int3d({$x#14739, $y#14738, $z#14737})), legion_domain_t(rect1d({$qBase#14740, (($qBase#14740+1000)-1)})))
      end
    end
  end
  var $particles_qSrcPart_25#14742 : partition#3844(disjoint, $particles_queue_25#14735, $primColors#14371) = partition(disjoint, $particles_queue_25#14735, $srcColoring#14736, $primColors#14371)
  legion_domain_point_coloring_destroy($srcColoring#14736)
  var $dstColoring#14743 : legion_domain_point_coloring_t = legion_domain_point_coloring_create()
  var $colorOff#14744 : int3d = int3d({-1, -1, -1})
  for $c#14745 : int3d($primColors#14371) in $primColors#14371 do
    var $srcBase#14746 : int64
    for $qptr#14747 : int1d(int8[376], $15627) in $particles_qSrcPart_25#14742[((($c#14745-$colorOff#14744)+{2, 3, 4})%{2, 3, 4})] do
      $srcBase#14746 = int64(int1d($qptr#14747))
      break
    end
    legion_domain_point_coloring_color_domain($dstColoring#14743, legion_domain_point_t($c#14745), legion_domain_t(rect1d({$srcBase#14746, (($srcBase#14746+1000)-1)})))
  end
  var $particles_qDstPart_25#14748 : partition#3846(aliased, $particles_queue_25#14735, $primColors#14371) = partition(aliased, $particles_queue_25#14735, $dstColoring#14743, $primColors#14371)
  legion_domain_point_coloring_destroy($dstColoring#14743)
  __parallelize_with $Fluid_primPart#14375, $particles_primPart#14383, $primColors#14371, (image($Fluid#14358, $particles_primPart#14383, $particles#14369.cell)<=$Fluid_primPart#14375) do
    particles_initValidField($particles#14369)
    Flow_InitializeCell($Fluid#14358, $Fluid#14358)
    Flow_InitializeCenterCoordinates($Fluid#14358, $Fluid#14358)
    Flow_InitializeUniform($Fluid#14358, $Fluid#14358)
    Flow_UpdateConservedFromPrimitive($Fluid#14358, $Fluid#14358)
    Flow_UpdateAuxiliaryVelocity($Fluid#14358, $Fluid#14358)
    Flow_UpdateGhostConservedStep1($Fluid#14358, $Fluid#14358)
    Flow_UpdateGhostConservedStep2($Fluid#14358, $Fluid#14358)
    Flow_UpdateGhostVelocityStep1($Fluid#14358, $Fluid#14358)
    Flow_UpdateGhostVelocityStep2($Fluid#14358, $Fluid#14358)
    Flow_ComputeVelocityGradientAll($Fluid#14358, $Fluid#14358)
    Flow_UpdateAuxiliaryThermodynamics($Fluid#14358, $Fluid#14358)
    Flow_UpdateGhostThermodynamicsStep1($Fluid#14358, $Fluid#14358)
    Flow_UpdateGhostThermodynamicsStep2($Fluid#14358, $Fluid#14358)
    Flow_UpdateGhostFieldsStep1($Fluid#14358, $Fluid#14358)
    Flow_UpdateGhostFieldsStep2($Fluid#14358, $Fluid#14358)
    $Particles_number#14356 = int64(int32(120))
    InitParticlesUniform($particles#14369, $Fluid#14358)
    $Flow_numberOfInteriorCells#14338 = int64(int32(0))
    $Flow_areaInterior#14341 = double(int32(0))
    $Flow_numberOfInteriorCells#14338 += numberOfInteriorCells($Fluid#14358, $Fluid#14358)
    $Flow_areaInterior#14341 += areaInterior($Fluid#14358, $Fluid#14358)
    $Flow_averagePressure#14352 = double(int32(0))
    $Flow_averageTemperature#14355 = double(int32(0))
    $Flow_averageKineticEnergy#14351 = double(int32(0))
    $Flow_minTemperature#14350 = double(int32(inf))
    $Flow_maxTemperature#14335 = double(int32(-inf))
    $Flow_averagePD#14343 = double(int32(0))
    $Flow_averageDissipation#14336 = double(int32(0))
    $Particles_averageTemperature#14342 = double(int32(0))
    $Flow_averagePressure#14352 += averagePressure($Fluid#14358, $Fluid#14358)
    $Flow_averageTemperature#14355 += averageTemperature($Fluid#14358, $Fluid#14358)
    $Flow_averageKineticEnergy#14351 += averageKineticEnergy($Fluid#14358, $Fluid#14358)
    $Flow_minTemperature#14350 min= minTemperature($Fluid#14358, $Fluid#14358)
    $Flow_maxTemperature#14335 max= maxTemperature($Fluid#14358, $Fluid#14358)
    $Particles_averageTemperature#14342 += Particles_IntegrateQuantities($particles#14369, $particles#14369)
    $Flow_averagePressure#14352 = ($Flow_averagePressure#14352/$Flow_areaInterior#14341)
    $Flow_averageTemperature#14355 = ($Flow_averageTemperature#14355/$Flow_areaInterior#14341)
    $Flow_averageKineticEnergy#14351 = ($Flow_averageKineticEnergy#14351/$Flow_areaInterior#14341)
    $Particles_averageTemperature#14342 = ($Particles_averageTemperature#14342/$Particles_number#14356)
    var $flag#14749 : int32 = int32((($TimeIntegrator_timeStep#14347%int32(333))==int32(0)))
    while ($flag#14749>0) do
      var $flag#14750 : int32 = int32((($TimeIntegrator_timeStep#14347%int32(555))==int32(0)))
      while ($flag#14750>0) do
        print($TimeIntegrator_deltaTime#14354)
        print_($Flow_minTemperature#14350, $Flow_maxTemperature#14335)
        print__($Particles_number#14356)
        print___()
        print____()
        $flag#14750 -= 1
      end
      print_____($TimeIntegrator_timeStep#14347, $TimeIntegrator_simTime#14353, $Flow_averagePressure#14352, $Flow_averageTemperature#14355, $Flow_averageKineticEnergy#14351, $Particles_averageTemperature#14342)
      $flag#14749 -= 1
    end
    var $flag#14751 : int32 = int32((($TimeIntegrator_timeStep#14347%int32(444))==int32(0)))
    while ($flag#14751>0) do
      var $filename#14752 : &int8 = &int8(malloc(uint64(256)))
      snprintf($filename#14752, uint64(256), "/home/manolis/proj/psaap/soleil-x/src/restart_fluid_%d.hdf", $TimeIntegrator_timeStep#14347)
      Fluid_hdf5create_rho_pressure_velocity($filename#14752)
      attach(hdf5, $Fluid_copy#14359.{rho, pressure, velocity}, $filename#14752, uint32(1))
      for $c#14753 : int3d($primColors#14371) in $primColors#14371 do
        var $p_r_c#14754 : region#3848(ispace#3790(int3d), Fluid_columns) = $Fluid_primPart#14375[int3d($c#14753)]
        var $p_s_c#14755 : region#3849(ispace#3791(int3d), Fluid_columns) = $Fluid_copy_primPart#14376[int3d($c#14753)]
        acquire($p_s_c#14755.{rho, pressure, velocity})
        copy($p_r_c#14754.{rho, pressure, velocity}, $p_s_c#14755.{rho, pressure, velocity})
        release($p_s_c#14755.{rho, pressure, velocity})
      end
      detach(hdf5, $Fluid_copy#14359.{rho, pressure, velocity})
      free(&opaque($filename#14752))
      $flag#14751 -= 1
    end
    var $flag#14756 : int32 = int32((($TimeIntegrator_timeStep#14347%int32(444))==int32(0)))
    while ($flag#14756>0) do
      var $filename#14757 : &int8 = &int8(malloc(uint64(256)))
      snprintf($filename#14757, uint64(256), "/home/manolis/proj/psaap/soleil-x/src/restart_particles_%d.hdf", $TimeIntegrator_timeStep#14347)
      particles_hdf5create_cell_position_particle_velocity_particle_temperature_diameter___valid($filename#14757)
      attach(hdf5, $particles_copy#14370.{cell, position, particle_velocity, particle_temperature, diameter, __valid}, $filename#14757, uint32(1))
      for $c#14758 : int3d($primColors#14371) in $primColors#14371 do
        var $p_r_c#14759 : region#3850(ispace#3792(int1d), particles_columns) = $particles_primPart#14383[int3d($c#14758)]
        var $p_s_c#14760 : region#3851(ispace#3793(int1d), particles_columns) = $particles_copy_primPart#14384[int3d($c#14758)]
        acquire($p_s_c#14760.{cell, position, particle_velocity, particle_temperature, diameter, __valid})
        copy($p_r_c#14759.{cell, position, particle_velocity, particle_temperature, diameter, __valid}, $p_s_c#14760.{cell, position, particle_velocity, particle_temperature, diameter, __valid})
        release($p_s_c#14760.{cell, position, particle_velocity, particle_temperature, diameter, __valid})
      end
      detach(hdf5, $particles_copy#14370.{cell, position, particle_velocity, particle_temperature, diameter, __valid})
      free(&opaque($filename#14757))
      $flag#14756 -= 1
    end
    while (($TimeIntegrator_simTime#14353<double(0.8029553)) and ($TimeIntegrator_timeStep#14347<int32(222))) do
      $maxC#14348 max= calculateConvectiveSpectralRadius($Fluid#14358, $Fluid#14358)
      $maxV#14345 max= calculateViscousSpectralRadius($Fluid#14358, $Fluid#14358)
      $maxH#14344 max= calculateHeatConductionSpectralRadius($Fluid#14358, $Fluid#14358)
      $TimeIntegrator_deltaTime#14354 = (double(0.4106583)/max($maxC#14348, max($maxV#14345, $maxH#14344)))
      Flow_InitializeTemporaries($Fluid#14358, $Fluid#14358)
      Particles_InitializeTemporaries($particles#14369, $particles#14369)
      $TimeIntegrator_timeOld#14337 = $TimeIntegrator_simTime#14353
      $TimeIntegrator_stage#14346 = int32(1)
      while ($TimeIntegrator_stage#14346<int32(5)) do
        Flow_InitializeTimeDerivatives($Fluid#14358, $Fluid#14358)
        Particles_InitializeTimeDerivatives($particles#14369, $particles#14369)
        Flow_UpdateGhostVelocityGradientStep1($Fluid#14358, $Fluid#14358)
        Flow_UpdateGhostVelocityGradientStep2($Fluid#14358, $Fluid#14358)
        Flow_AddGetFlux($Fluid#14358, $Fluid#14358)
        Flow_AddUpdateUsingFlux($Fluid#14358, $Fluid#14358)
        Flow_AddBodyForces($Fluid#14358, $Fluid#14358)
        Particles_LocateInCells($particles#14369, $particles#14369)
        for $c#14761 : int3d($primColors#14371) in $primColors#14371 do
          particles_pushAll(int3d($c#14761), $particles_primPart#14383[int3d($c#14761)], $particles_qSrcPart_0#14392[int3d($c#14761)], $particles_qSrcPart_1#14406[int3d($c#14761)], $particles_qSrcPart_2#14420[int3d($c#14761)], $particles_qSrcPart_3#14434[int3d($c#14761)], $particles_qSrcPart_4#14448[int3d($c#14761)], $particles_qSrcPart_5#14462[int3d($c#14761)], $particles_qSrcPart_6#14476[int3d($c#14761)], $particles_qSrcPart_7#14490[int3d($c#14761)], $particles_qSrcPart_8#14504[int3d($c#14761)], $particles_qSrcPart_9#14518[int3d($c#14761)], $particles_qSrcPart_10#14532[int3d($c#14761)], $particles_qSrcPart_11#14546[int3d($c#14761)], $particles_qSrcPart_12#14560[int3d($c#14761)], $particles_qSrcPart_13#14574[int3d($c#14761)], $particles_qSrcPart_14#14588[int3d($c#14761)], $particles_qSrcPart_15#14602[int3d($c#14761)], $particles_qSrcPart_16#14616[int3d($c#14761)], $particles_qSrcPart_17#14630[int3d($c#14761)], $particles_qSrcPart_18#14644[int3d($c#14761)], $particles_qSrcPart_19#14658[int3d($c#14761)], $particles_qSrcPart_20#14672[int3d($c#14761)], $particles_qSrcPart_21#14686[int3d($c#14761)], $particles_qSrcPart_22#14700[int3d($c#14761)], $particles_qSrcPart_23#14714[int3d($c#14761)], $particles_qSrcPart_24#14728[int3d($c#14761)], $particles_qSrcPart_25#14742[int3d($c#14761)])
        end
        for $c#14762 : int3d($primColors#14371) in $primColors#14371 do
          particles_pullAll(int3d($c#14762), $particles_primPart#14383[int3d($c#14762)], $particles_qDstPart_0#14398[int3d($c#14762)], $particles_qDstPart_1#14412[int3d($c#14762)], $particles_qDstPart_2#14426[int3d($c#14762)], $particles_qDstPart_3#14440[int3d($c#14762)], $particles_qDstPart_4#14454[int3d($c#14762)], $particles_qDstPart_5#14468[int3d($c#14762)], $particles_qDstPart_6#14482[int3d($c#14762)], $particles_qDstPart_7#14496[int3d($c#14762)], $particles_qDstPart_8#14510[int3d($c#14762)], $particles_qDstPart_9#14524[int3d($c#14762)], $particles_qDstPart_10#14538[int3d($c#14762)], $particles_qDstPart_11#14552[int3d($c#14762)], $particles_qDstPart_12#14566[int3d($c#14762)], $particles_qDstPart_13#14580[int3d($c#14762)], $particles_qDstPart_14#14594[int3d($c#14762)], $particles_qDstPart_15#14608[int3d($c#14762)], $particles_qDstPart_16#14622[int3d($c#14762)], $particles_qDstPart_17#14636[int3d($c#14762)], $particles_qDstPart_18#14650[int3d($c#14762)], $particles_qDstPart_19#14664[int3d($c#14762)], $particles_qDstPart_20#14678[int3d($c#14762)], $particles_qDstPart_21#14692[int3d($c#14762)], $particles_qDstPart_22#14706[int3d($c#14762)], $particles_qDstPart_23#14720[int3d($c#14762)], $particles_qDstPart_24#14734[int3d($c#14762)], $particles_qDstPart_25#14748[int3d($c#14762)])
        end
        Particles_AddFlowCoupling($particles#14369, $particles#14369, $Fluid#14358)
        Particles_AddBodyForces($particles#14369, $particles#14369)
        AddRadiation($particles#14369)
        $Flow_averageHeatSource#14334 += Flow_AddParticlesCoupling($particles#14369, $particles#14369, $Fluid#14358)
        Flow_UpdateVars($Fluid#14358, $Fluid#14358, $TimeIntegrator_deltaTime#14354, $TimeIntegrator_stage#14346)
        Particles_UpdateVars($particles#14369, $particles#14369, $TimeIntegrator_deltaTime#14354, $TimeIntegrator_stage#14346)
        Flow_UpdateAuxiliaryVelocity($Fluid#14358, $Fluid#14358)
        Flow_UpdateGhostConservedStep1($Fluid#14358, $Fluid#14358)
        Flow_UpdateGhostConservedStep2($Fluid#14358, $Fluid#14358)
        Flow_UpdateGhostVelocityStep1($Fluid#14358, $Fluid#14358)
        Flow_UpdateGhostVelocityStep2($Fluid#14358, $Fluid#14358)
        Flow_ComputeVelocityGradientAll($Fluid#14358, $Fluid#14358)
        Flow_UpdateAuxiliaryThermodynamics($Fluid#14358, $Fluid#14358)
        Flow_UpdateGhostThermodynamicsStep1($Fluid#14358, $Fluid#14358)
        Flow_UpdateGhostThermodynamicsStep2($Fluid#14358, $Fluid#14358)
        Particles_UpdateAuxiliaryStep1($particles#14369, $particles#14369)
        Particles_UpdateAuxiliaryStep2($particles#14369, $particles#14369)
        $TimeIntegrator_simTime#14353 = ($TimeIntegrator_timeOld#14337+((double(0.5)*(int32(1)+($TimeIntegrator_stage#14346/int32(3))))*$TimeIntegrator_deltaTime#14354))
        $TimeIntegrator_stage#14346 = ($TimeIntegrator_stage#14346+int32(1))
        for $c#14763 : int3d($primColors#14371) in $primColors#14371 do
          $Particles_number#14356 += Particles_DeleteEscapingParticles($particles_primPart#14383[int3d($c#14763)], $particles_primPart#14383[int3d($c#14763)])
        end
      end
      $TimeIntegrator_timeStep#14347 = ($TimeIntegrator_timeStep#14347+int32(1))
      var $flag#14764 : int32 = int32((($TimeIntegrator_timeStep#14347%int32(333))==int32(0)))
      while ($flag#14764>0) do
        $Flow_averagePressure#14352 = double(int32(0))
        $Flow_averageTemperature#14355 = double(int32(0))
        $Flow_averageKineticEnergy#14351 = double(int32(0))
        $Flow_minTemperature#14350 = double(int32(inf))
        $Flow_maxTemperature#14335 = double(int32(-inf))
        $Flow_averagePD#14343 = double(int32(0))
        $Flow_averageDissipation#14336 = double(int32(0))
        $Particles_averageTemperature#14342 = double(int32(0))
        $Flow_averagePressure#14352 += averagePressure($Fluid#14358, $Fluid#14358)
        $Flow_averageTemperature#14355 += averageTemperature($Fluid#14358, $Fluid#14358)
        $Flow_averageKineticEnergy#14351 += averageKineticEnergy($Fluid#14358, $Fluid#14358)
        $Flow_minTemperature#14350 min= minTemperature($Fluid#14358, $Fluid#14358)
        $Flow_maxTemperature#14335 max= maxTemperature($Fluid#14358, $Fluid#14358)
        $Particles_averageTemperature#14342 += Particles_IntegrateQuantities($particles#14369, $particles#14369)
        $Flow_averagePressure#14352 = ($Flow_averagePressure#14352/$Flow_areaInterior#14341)
        $Flow_averageTemperature#14355 = ($Flow_averageTemperature#14355/$Flow_areaInterior#14341)
        $Flow_averageKineticEnergy#14351 = ($Flow_averageKineticEnergy#14351/$Flow_areaInterior#14341)
        $Particles_averageTemperature#14342 = ($Particles_averageTemperature#14342/$Particles_number#14356)
        var $flag#14765 : int32 = int32((($TimeIntegrator_timeStep#14347%int32(333))==int32(0)))
        while ($flag#14765>0) do
          var $flag#14766 : int32 = int32((($TimeIntegrator_timeStep#14347%int32(555))==int32(0)))
          while ($flag#14766>0) do
            print______($TimeIntegrator_deltaTime#14354)
            print_______($Flow_minTemperature#14350, $Flow_maxTemperature#14335)
            print________($Particles_number#14356)
            print_________()
            print__________()
            $flag#14766 -= 1
          end
          print___________($TimeIntegrator_timeStep#14347, $TimeIntegrator_simTime#14353, $Flow_averagePressure#14352, $Flow_averageTemperature#14355, $Flow_averageKineticEnergy#14351, $Particles_averageTemperature#14342)
          $flag#14765 -= 1
        end
        var $flag#14767 : int32 = int32((($TimeIntegrator_timeStep#14347%int32(444))==int32(0)))
        while ($flag#14767>0) do
          var $filename#14768 : &int8 = &int8(malloc(uint64(256)))
          snprintf($filename#14768, uint64(256), "/home/manolis/proj/psaap/soleil-x/src/restart_fluid_%d.hdf", $TimeIntegrator_timeStep#14347)
          Fluid_hdf5create_rho_pressure_velocity_($filename#14768)
          attach(hdf5, $Fluid_copy#14359.{rho, pressure, velocity}, $filename#14768, uint32(1))
          for $c#14769 : int3d($primColors#14371) in $primColors#14371 do
            var $p_r_c#14770 : region#3907(ispace#3849(int3d), Fluid_columns) = $Fluid_primPart#14375[int3d($c#14769)]
            var $p_s_c#14771 : region#3908(ispace#3850(int3d), Fluid_columns) = $Fluid_copy_primPart#14376[int3d($c#14769)]
            acquire($p_s_c#14771.{rho, pressure, velocity})
            copy($p_r_c#14770.{rho, pressure, velocity}, $p_s_c#14771.{rho, pressure, velocity})
            release($p_s_c#14771.{rho, pressure, velocity})
          end
          detach(hdf5, $Fluid_copy#14359.{rho, pressure, velocity})
          free(&opaque($filename#14768))
          $flag#14767 -= 1
        end
        var $flag#14772 : int32 = int32((($TimeIntegrator_timeStep#14347%int32(444))==int32(0)))
        while ($flag#14772>0) do
          var $filename#14773 : &int8 = &int8(malloc(uint64(256)))
          snprintf($filename#14773, uint64(256), "/home/manolis/proj/psaap/soleil-x/src/restart_particles_%d.hdf", $TimeIntegrator_timeStep#14347)
          particles_hdf5create_cell_position_particle_velocity_particle_temperature_diameter___valid_($filename#14773)
          attach(hdf5, $particles_copy#14370.{cell, position, particle_velocity, particle_temperature, diameter, __valid}, $filename#14773, uint32(1))
          for $c#14774 : int3d($primColors#14371) in $primColors#14371 do
            var $p_r_c#14775 : region#3909(ispace#3851(int1d), particles_columns) = $particles_primPart#14383[int3d($c#14774)]
            var $p_s_c#14776 : region#3910(ispace#3852(int1d), particles_columns) = $particles_copy_primPart#14384[int3d($c#14774)]
            acquire($p_s_c#14776.{cell, position, particle_velocity, particle_temperature, diameter, __valid})
            copy($p_r_c#14775.{cell, position, particle_velocity, particle_temperature, diameter, __valid}, $p_s_c#14776.{cell, position, particle_velocity, particle_temperature, diameter, __valid})
            release($p_s_c#14776.{cell, position, particle_velocity, particle_temperature, diameter, __valid})
          end
          detach(hdf5, $particles_copy#14370.{cell, position, particle_velocity, particle_temperature, diameter, __valid})
          free(&opaque($filename#14773))
          $flag#14772 -= 1
        end
        $flag#14764 -= 1
      end
    end
    var $flag#14777 : int32 = int32((($TimeIntegrator_timeStep#14347%int32(333))==int32(0)))
    while ($flag#14777>0) do
      var $flag#14778 : int32 = int32((($TimeIntegrator_timeStep#14347%int32(555))==int32(0)))
      while ($flag#14778>0) do
        print____________($TimeIntegrator_deltaTime#14354)
        print_____________($Flow_minTemperature#14350, $Flow_maxTemperature#14335)
        print______________($Particles_number#14356)
        print_______________()
        print________________()
        $flag#14778 -= 1
      end
      print_________________($TimeIntegrator_timeStep#14347, $TimeIntegrator_simTime#14353, $Flow_averagePressure#14352, $Flow_averageTemperature#14355, $Flow_averageKineticEnergy#14351, $Particles_averageTemperature#14342)
      $flag#14777 -= 1
    end
  end
end
