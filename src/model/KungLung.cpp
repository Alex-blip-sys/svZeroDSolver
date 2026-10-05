// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "KungLung.h"

void KungLung::setup_dofs(DOFHandler& dofhandler) {
  Block::setup_dofs_(
      dofhandler, 5,
      {"pressure_top", "pressure_bottom", "pressure_mid"});
}

void KungLung::update_constant(SparseSystem& system,
                               std::vector<double>& parameters) {
  const double Rpa_top =
      parameters[global_param_ids[ParamId::RPA_TOP]];
  const double Rpa_bottom =
      parameters[global_param_ids[ParamId::RPA_BOTTOM]];
  const double Rpv_top =
      parameters[global_param_ids[ParamId::RPV_TOP]];
  const double Rpv_bottom =
      parameters[global_param_ids[ParamId::RPV_BOTTOM]];
  const double C_top =
      parameters[global_param_ids[ParamId::C_TOP]];
  const double C_bottom =
      parameters[global_param_ids[ParamId::C_BOTTOM]];

  // Local variable ordering:
  // 0 -> P_in
  // 1 -> Q_in
  // 2 -> P_out
  // 3 -> Q_out
  // 4 -> P_top
  // 5 -> P_bottom
  // 6 -> P_mid

  // ------------------------------------------------------------
  // Eq 0:
  // (P_in - P_top)/Rpa_top
  // - (P_top - P_out)/Rpv_top
  // - C_top * (dP_top/dt - dP_mid/dt) = 0
  // ------------------------------------------------------------
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0 / Rpa_top;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[2]) = 1.0 / Rpv_top;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[4]) =
      -(1.0 / Rpa_top + 1.0 / Rpv_top);

  system.E.coeffRef(global_eqn_ids[0], global_var_ids[4]) = -C_top;
  system.E.coeffRef(global_eqn_ids[0], global_var_ids[6]) =  C_top;

  // ------------------------------------------------------------
  // Eq 1:
  // (P_in - P_bottom)/Rpa_bottom
  // - (P_bottom - P_out)/Rpv_bottom
  // - C_bottom * (dP_bottom/dt - dP_mid/dt) = 0
  // ------------------------------------------------------------
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[0]) = 1.0 / Rpa_bottom;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[2]) = 1.0 / Rpv_bottom;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[5]) =
      -(1.0 / Rpa_bottom + 1.0 / Rpv_bottom);

  system.E.coeffRef(global_eqn_ids[1], global_var_ids[5]) = -C_bottom;
  system.E.coeffRef(global_eqn_ids[1], global_var_ids[6]) =  C_bottom;

  // ------------------------------------------------------------
  // Eq 2:
  // Q_in - (P_in - P_top)/Rpa_top - (P_in - P_bottom)/Rpa_bottom = 0
  // ------------------------------------------------------------
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[0]) =
      -(1.0 / Rpa_top + 1.0 / Rpa_bottom);
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[1]) = 1.0;
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[4]) = 1.0 / Rpa_top;
  system.F.coeffRef(global_eqn_ids[2], global_var_ids[5]) = 1.0 / Rpa_bottom;

  // ------------------------------------------------------------
  // Eq 3:
  // Q_out - Q_in = 0
  // ------------------------------------------------------------
  system.F.coeffRef(global_eqn_ids[3], global_var_ids[1]) = -1.0;
  system.F.coeffRef(global_eqn_ids[3], global_var_ids[3]) =  1.0;

  // ------------------------------------------------------------
  // Eq 4:
  // P_mid - P_out - P_source(t) = 0
  // ------------------------------------------------------------
  system.F.coeffRef(global_eqn_ids[4], global_var_ids[2]) = -1.0;
  system.F.coeffRef(global_eqn_ids[4], global_var_ids[6]) =  1.0;
}

void KungLung::update_time(SparseSystem& system,
                           std::vector<double>& parameters) {
  const double Psrc =
      parameters[global_param_ids[ParamId::P_SOURCE]];

  // IMPORTANT:
  // Use the same style here that your local PressureReferenceBC::update_time()
  // uses for the forcing vector. In many builds this is:
  system.C(global_eqn_ids[4]) = -Psrc;

  // If your SparseSystem uses a different member name for the forcing vector,
  // mirror PressureReferenceBC::update_time() exactly and assign -Psrc to Eq 4.
}
