/*
 * This file is part of 1Dchain.
 *
 * Copyright (C) 2026
 * Vincenzo Di Florio, Davide Carbone
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#ifndef PARAM_H
#define PARAM_H

#include <unordered_map>
#include <iostream>
#include <string>

struct
  params_t {

  std::unordered_map<std::string, int> iparams;
  std::unordered_map<std::string, double> dparams;
  std::unordered_map<std::string, std::string> sparams;
  std::unordered_map<std::string, bool> bparams;

  void
  read (const std::string &fname);

  void
  write (const std::string &fname);
};

std::ostream&
operator<< (std::ostream& os, const params_t& p);

std::istream&
operator>> (std::istream& is, params_t& p);

extern params_t p;

#endif
