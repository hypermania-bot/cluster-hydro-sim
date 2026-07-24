#ifndef IO_HPP
#define IO_HPP
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <optional>



template<typename Scalar>
void write_to_file(const std::vector<Scalar> &vector, const std::string &filename)
{
  //char *memblock = (char *)&vector[0];
  const char *memblock = reinterpret_cast<const char *>(vector.data());
  std::ofstream file(filename, std::ios::binary);
  if(file.is_open()){
    file.write(memblock, vector.size() * sizeof(Scalar));
  }
}

#ifndef NOT_USING_EIGEN

#include <Eigen/Dense>
template<typename Derived>
void write_to_file(const Eigen::PlainObjectBase<Derived> &obj, const std::string &filename)
{
  std::ofstream file(filename, std::ios::binary);
  if(file.is_open()){
    file.write((char *)obj.data(), obj.size() * sizeof(typename Eigen::PlainObjectBase<Derived>::Scalar));
  }
}

template<typename Derived>
void write_to_filename_template(const Eigen::PlainObjectBase<Derived> &obj, const std::string &format_string, const int idx)
{
  char filename[128];
  sprintf(filename, format_string.data(), idx);
  std::ofstream file(filename, std::ios::binary);
  if(file.is_open()){
    file.write((char *)obj.data(), obj.size() * sizeof(typename Eigen::PlainObjectBase<Derived>::Scalar));
  }
}

#endif

template<typename ReturnType>
std::optional<ReturnType> read_from_file(const std::string &filename) {
  std::ifstream file(filename, std::ios::in | std::ios::binary | std::ios::ate);
  if(file.is_open()) {
    const std::streampos size = file.tellg();
    file.seekg(0, std::ios::beg);
    
    // Assume ReturnType consists of double by default
    constexpr size_t elem_size = sizeof(decltype(ReturnType()[0]));
    // std::cout << "(read_from_file) deduced elem_size = " << elem_size << std::endl;
    const size_t obj_size = (static_cast<size_t>(size) + elem_size - 1) / elem_size;
    std::optional<ReturnType> result{ReturnType(obj_size)};
    file.read(reinterpret_cast<char *>(result->data()), size);
    return result;
  } else {
    return std::optional<ReturnType>();
  }
}

template<typename ReturnType>
std::optional<ReturnType> read_from_file_template(const std::string &format_string, const int idx) {
  char filename[128];
  sprintf(filename, format_string.data(), idx);
  return read_from_file<ReturnType>(filename);
}


#endif
