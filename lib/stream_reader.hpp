#ifndef GLVIS_STREAM_READER_HPP
#define GLVIS_STREAM_READER_HPP

#include <string>
#include "mfem.hpp"

struct StreamState
{
   mfem::Vector sol, solu, solv, solw, normals;
   std::string keys;
   mfem::Mesh *mesh{nullptr};
   mfem::GridFunction *grid_f{nullptr};
   int is_gf{0};
   bool fix_elem_orient{false};
   bool save_coloring{false};
};
extern StreamState stream_state;

void Extrude1DMeshAndSolution(mfem::Mesh **mesh_p,
                              mfem::GridFunction **grid_f_p,
                              mfem::Vector *sol);
void SetMeshSolution(mfem::Mesh *mesh, mfem::GridFunction *&grid_f,
                     bool save_coloring);
int ReadStream(std::istream &is, const std::string &data_type);

#endif // GLVIS_STREAM_READER_HPP
