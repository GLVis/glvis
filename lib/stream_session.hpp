// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-443271.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see https://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.
#pragma once

#include <memory>
#include <string>
#include <vector>

int GLVisStreamSession(const bool fix_elem_orient,
                       const bool save_coloring,
                       const bool keep_attr,
                       const bool headless,
                       const std::string &plot_caption,
                       const std::string &data_type,
                       std::vector<std::unique_ptr<std::istream>> &&streams);