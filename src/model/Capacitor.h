// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause
#ifndef SVZERODSOLVER_MODEL_CAPACITOR_HPP_
#define SVZERODSOLVER_MODEL_CAPACITOR_HPP_

#include "Block.h"
#include "SparseSystem.h"

/**
 * @brief Flow-through RC windkessel for pulmonary circulation.
 *
 * Models a capacitor where the inlet and outlet
 * flows are equal (Q_in = Q_out). The capacitor stores pressure but does
 * not store blood volume. This matches the pulmonary model in
 * \cite sankaran2012patient (equation X6' = X4 - X7).
 *
 * ### Governing equations
 *
 * \f[
 * C_p \, (\dot{P}_{in} - \dot{P}_{out}) - Q_{in} = 0
 * \f]
 * \f[
 * Q_{in} - Q_{out} = 0
 * \f]
 *
 * ### Parameters
 *
 * * `Cp` — Pulmonary capacitance
 */
class Capacitor : public Block {
 public:
  /**
   * @brief Construct a new Capacitor object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  Capacitor(int id, Model* model)
      : Block(id, model, BlockType::capacitor, BlockClass::vessel,
              {{"Cp", InputParameter()}}) {}

  /// @brief Local IDs of the parameters
  enum ParamId { CP = 0 };

  void setup_dofs(DOFHandler& dofhandler) override;
  void update_constant(SparseSystem& system,
                       std::vector<double>& parameters) override;

  /// @brief Number of triplets of element
  TripletsContributions num_triplets{3, 2, 0};
};

#endif  // SVZERODSOLVER_MODEL_CAPACITOR_HPP_
