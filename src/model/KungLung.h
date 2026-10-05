// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#ifndef SVZERODSOLVER_MODEL_KUNGLUNG_HPP_
#define SVZERODSOLVER_MODEL_KUNGLUNG_HPP_

#include "Block.h"
#include "SparseSystem.h"

/**
 * @brief Lung element based on the exercise Fontan pulmonary circuit.
 *
 * This block represents two parallel arterial resistors, two parallel venous
 * resistors, two capacitors coupled to an internal midpoint, and a prescribed
 * time-varying pressure source applied between the capacitor midpoint and the
 * outlet node.
 *
 * External variables:
 *   P_in, Q_in, P_out, Q_out
 *
 * Internal variables:
 *   pressure_top    = P_top
 *   pressure_bottom = P_bot
 *   pressure_mid    = P_mid
 *
 * The source waveform is prescribed through:
 *   t_source, P_source
 */
class KungLung : public Block {
 public:
  enum ParamId {
    RPA_TOP = 0,
    RPA_BOTTOM = 1,
    RPV_TOP = 2,
    RPV_BOTTOM = 3,
    C_TOP = 4,
    C_BOTTOM = 5,
    T_SOURCE = 6,
    P_SOURCE = 7
  };

  KungLung(int id, Model* model)
      : Block(id, model, BlockType::kung_lung, BlockClass::vessel,
              {
                  {"Rpa_top", InputParameter()},
                  {"Rpa_bottom", InputParameter()},
                  {"Rpv_top", InputParameter()},
                  {"Rpv_bottom", InputParameter()},
                  {"C_top", InputParameter()},
                  {"C_bottom", InputParameter()},
                  {"t_source", InputParameter(false, true)},
                  {"P_source", InputParameter(false, true)},
              }) {}

  void setup_dofs(DOFHandler& dofhandler) override;

  void update_constant(SparseSystem& system,
                       std::vector<double>& parameters) override;

  void update_time(SparseSystem& system,
                   std::vector<double>& parameters) override;

  /// F entries = 14, E entries = 4, D entries = 0
  TripletsContributions num_triplets{14, 4, 0};
};

#endif  // SVZERODSOLVER_MODEL_KUNGLUNG_HPP_
