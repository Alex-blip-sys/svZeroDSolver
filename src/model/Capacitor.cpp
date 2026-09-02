// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#include "Capacitor.h"

void Capacitor::setup_dofs(DOFHandler& dofhandler) {
  Block::setup_dofs_(dofhandler, 1, {});
}

void Capacitor::update_constant(SparseSystem& system,
                                    std::vector<double>& parameters) {
  double Rpd = parameters[global_param_ids[ParamId::RPD]];
  double Cp = parameters[global_param_ids[ParamId::CP]];

  // Eq 0: Cp*dP_in - Q_in = 0
  system.E.coeffRef(global_eqn_ids[0], global_var_ids[0]) = Cp;
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) = -1.0;

  // Eq 1: Q_in - Q_out = 0
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[1]) = 1.0;
  system.F.coeffRef(global_eqn_ids[1], global_var_ids[3]) = -1.0;
}
