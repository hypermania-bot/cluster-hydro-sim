/*!
  \file param.hpp
  \brief Utilities for managing parameters specifying different kinds of simulations.
  
  Data types and utilities for managing the parameters of a simulation.
  The parameters are collected as a struct of double or long long int.
  Uses Boost.PFR for saving struct member type names, 
  so that Mathematica can import the struct in binary.
*/
#ifndef PARAM_HPP
#define PARAM_HPP

#include "utility.hpp"
#include "boost/pfr.hpp"
#include "boost/type_index.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <boost/multiprecision/float128.hpp>


struct SampleParam {
  long long int N;
  double L;
  double m;
  double lambda;
  double k_ast;
  double k_Psi;
  double varphi_std_dev;
  double Psi_std_dev;
  double a1;
  double H1;
  double t1;
};


template<typename T>
void print_param(const T &param) {
  auto names = boost::pfr::names_as_array<T>();
  auto func = [&](const auto &field, std::size_t i) {
		std::cout << names[i] << ": " << field
			  << " (" << boost::typeindex::type_id_runtime(field) << ")\n";
	      };

  auto c = [&](){ boost::pfr::for_each_field(param, func); };
  run_and_print("The parameters for the simulation", c);
}


template<typename T>
void save_param_names(const std::string &filename) {
  std::ofstream outstream(filename);
  auto names = boost::pfr::names_as_array<T>();
  for(auto name : names) {
    outstream << name << '\n';
  }
}

/*
// Compiles with Intel icpx, but doesn't compile with gcc due to "Explicit template specialization cannot have a storage class"
template<typename T> std::string_view Mathematica_format;

template<> constexpr static std::string_view Mathematica_format<double> = "Real64";

template<> constexpr static std::string_view Mathematica_format<long long int> = "Integer64";
*/

/*
// Compiles with gcc, fails at link stage with Intel icpx due to multiple definitions
template<typename T> std::string_view Mathematica_format;

template<> constexpr std::string_view Mathematica_format<double> = "Real64";

template<> constexpr std::string_view Mathematica_format<long long int> = "Integer64";
*/

/*
// Alternate solution
constexpr static std::string get_Mathematica_format(double) {
  return std::string("Real64");
}

constexpr static std::string get_Mathematica_format(long long int) {
  return std::string("Integer64");
}
*/

namespace {
  template<typename T> std::string_view Mathematica_format;
  template<> constexpr std::string_view Mathematica_format<double> = "Real64";
  template<> constexpr std::string_view Mathematica_format<boost::multiprecision::float128> = "Real128";
  template<> constexpr std::string_view Mathematica_format<long long int> = "Integer64";
  template<size_t size> constexpr std::string_view Mathematica_format<std::array<double, size>> = "Real64";
}

template<typename T>
void save_param_Mathematica_formats(const std::string &filename) {
  std::ofstream outstream(filename);
  auto func = [&](const auto &field) {
		typedef std::remove_const_t<std::remove_reference_t<decltype(field)>> type_of_field;
		outstream << Mathematica_format<type_of_field> << '\n';
	      };
  boost::pfr::for_each_field(T(), func);
}


template<typename T>
static void save_param(const T &param, const std::string &filename){
  std::ofstream outstream(filename, std::ios::binary);
  if(outstream.is_open()){
    outstream.write((const char *)&param, sizeof(T));
  }
}


template<typename T>
void save_param_offsets(const std::string &filename) {
  std::vector<size_t> offsets;
  T t;
  T *t_ptr = &t;
  auto func = [&](const auto &field) {
    auto *field_ptr = &field;
    offsets.push_back(reinterpret_cast<size_t>(field_ptr) - reinterpret_cast<size_t>(t_ptr));
  };
  boost::pfr::for_each_field(t, func);
  
  std::ofstream outstream(filename);  
  for(auto offset : offsets) {
    outstream << offset << '\n';
  }
}


template<typename T>
void save_param_for_Mathematica(const T &param, const std::string &dir) {
  save_param_names<T>(dir + "paramNames.txt");
  save_param_Mathematica_formats<T>(dir + "paramTypes.txt");
  save_param<T>(param, dir + "param.dat");
  save_param_offsets<T>(dir + "paramOffsets.txt");
}


template<typename T>
void save_param_types(const std::string &filename) {
  std::ofstream outstream(filename);
  auto func = [&](const auto &field) {
		outstream << boost::typeindex::type_id_runtime(field) << '\n';
	      };
  boost::pfr::for_each_field(T(), func);
}






#endif
