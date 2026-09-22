(* Export structural nonzeros without evaluating at a special physical state.
   Run moving_three_fluid_compiler.wl first to regenerate its coefficients. *)
directory = DirectoryName[$InputFileName];
coefficients = Get[FileNameJoin[{directory, "moving_three_fluid_coefficients.wl"}]];
Export[FileNameJoin[{directory, "moving_band_pattern.json"}],
 AssociationMap[Map[If[# === 0, 0, 1] &,
     coefficients[#]["Matrix"], {2}] &, Keys[coefficients]], "RawJSON"];
