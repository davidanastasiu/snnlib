/**
 * Copyright (c) David C. Anastasiu
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */

// -*- c++ -*-

#pragma once

#include <cstdint>

#define SNNLIB_VERSION_MAJOR 0
#define SNNLIB_VERSION_MINOR 1
#define SNNLIB_VERSION_PATCH 0

namespace snnlib
{
    static constexpr struct {
        uint8_t major, minor, patch;
    } version = {
            SNNLIB_VERSION_MAJOR,
            SNNLIB_VERSION_MINOR,
            SNNLIB_VERSION_PATCH
    };
} /* end namespace */
