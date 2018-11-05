#include "visual.hpp"
#include "palettes.hpp"
#include <SDL2/SDL_hints.h>
#include <emscripten/bind.h>

std::string plot_caption;
std::string extra_caption;
mfem::GeometryRefiner GLVisGeometryRefiner;

namespace js
{

using namespace mfem;

struct vis_data {
    int field_type;
    mfem::Mesh * mesh;
    mfem::Vector sol,
        solu,
        solv,
        solw,
        normals;
    mfem::GridFunction * grid_f;
    std::string keys;
} glvis_data;

void Extrude1DMeshAndSolution(vis_data & data) {
    if (data.mesh->Dimension() != 1 || data.mesh->SpaceDimension() != 1)
    {
        return;
    }

    // find xmin and xmax over the vertices of the 1D mesh
    double xmin = numeric_limits<double>::infinity();
    double xmax = -xmin;
    for (int i = 0; i < data.mesh->GetNV(); i++)
    {
        const double x = data.mesh->GetVertex(i)[0];
        if (x < xmin)
        {
            xmin = x;
        }
        if (x > xmax)
        {
            xmax = x;
        }
    }
    Mesh *mesh2d = Extrude1D(data.mesh, 1, 0.1*(xmax - xmin));

    if (data.grid_f)
    {
        GridFunction *grid_f_2d =
            Extrude1DGridFunction(data.mesh, mesh2d, data.grid_f, 1);
        delete data.grid_f;
        data.grid_f = grid_f_2d;
    }
    else if (data.sol.Size() == data.mesh->GetNV())
    {
        Vector sol2d(mesh2d->GetNV());
        for (int i = 0; i < data.mesh->GetNV(); i++)
        {
            sol2d(2*i+0) = sol2d(2*i+1) = data.sol(i);
        }
        data.sol = sol2d;
    }
    data.mesh = mesh2d;
}

void setMeshSolution(vis_data & data)
{
   if (1) // checkerboard solution
   {
      FiniteElementCollection *cfec;
      if (data.mesh->Dimension() == 1)
      {
         cfec = new L2_FECollection(0, 1);
      }
      else if (data.mesh->Dimension() == 2)
      {
         cfec = new Const2DFECollection;
      }
      else
      {
         cfec = new Const3DFECollection;
      }
      FiniteElementSpace *cfes = new FiniteElementSpace(data.mesh, cfec);
      data.grid_f = new GridFunction(cfes);
      data.grid_f->MakeOwner(cfec);
      {
         Array<int> coloring;
         srandom(time(0));
         double a = double(random()) / (double(RAND_MAX) + 1.);
         int el0 = (int)floor(a * data.mesh->GetNE());
         cout << "Generating coloring starting with element " << el0+1
              << " / " << data.mesh->GetNE() << endl;
         data.mesh->GetElementColoring(coloring, el0);
         for (int i = 0; i < coloring.Size(); i++)
         {
            (*data.grid_f)(i) = coloring[i];
         }
         cout << "Number of colors: " << data.grid_f->Max() + 1 << endl;
      }
      data.grid_f->GetNodalValues(data.sol);
      //is_gf = 1;
   }
   else // zero solution
   {
      data.sol.SetSize (data.mesh->GetNV());
      data.sol = 0.0;
   }
}

vis_data readStream(std::string& input, std::string& data_type) {
    std::stringstream mesh_file;
    mesh_file.str(input);
    mesh_file.clear();
    bool fix_elem_orient = false;
    vis_data data;
    data.mesh = nullptr;
    data.grid_f = nullptr;
    if (data_type == "fem2d_data")
    {
        data.mesh = new mfem::Mesh(mesh_file, 0, 0, fix_elem_orient);
        data.sol.Load(mesh_file, data.mesh->GetNV());
    }
    else if (data_type == "vfem2d_data" || data_type == "vfem2d_data_keys")
    {
        data.field_type = 1;
        data.mesh = new mfem::Mesh(mesh_file, 0, 0, fix_elem_orient);
        data.solu.Load(mesh_file, data.mesh->GetNV());
        data.solv.Load(mesh_file, data.mesh->GetNV());
        if (data_type == "vfem2d_data_keys")
        {
            mesh_file >> data.keys;
        }
    }
    else if (data_type == "fem3d_data")
    {
        data.mesh = new mfem::Mesh(mesh_file, 0, 0, fix_elem_orient);
        data.sol.Load(mesh_file, data.mesh->GetNV());
    }
    else if (data_type == "vfem3d_data" || data_type == "vfem3d_data_keys")
    {
        data.field_type = 1;
        data.mesh = new mfem::Mesh(mesh_file, 0, 0, fix_elem_orient);
        data.solu.Load(mesh_file, data.mesh->GetNV());
        data.solv.Load(mesh_file, data.mesh->GetNV());
        data.solw.Load(mesh_file, data.mesh->GetNV());
        if (data_type == "vfem3d_data_keys")
        {
            mesh_file >> data.keys;
        }
    }
    else if (data_type == "fem2d_gf_data" || data_type == "fem2d_gf_data_keys")
    {
        data.mesh = new mfem::Mesh(mesh_file, 1, 0, fix_elem_orient);
        data.grid_f = new mfem::GridFunction(data.mesh, mesh_file);
        if (data_type == "fem2d_gf_data_keys")
        {
            mesh_file >> data.keys;
        }
    }
    else if (data_type == "vfem2d_gf_data" || data_type == "vfem2d_gf_data_keys")
    {
        data.field_type = 1;
        data.mesh = new mfem::Mesh(mesh_file, 1, 0, fix_elem_orient);
        data.grid_f = new mfem::GridFunction(data.mesh, mesh_file);
        if (data_type == "vfem2d_gf_data_keys")
        {
            mesh_file >> data.keys;
        }
    }
    else if (data_type == "fem3d_gf_data" || data_type == "fem3d_gf_data_keys")
    {
        data.mesh = new mfem::Mesh(mesh_file, 1, 0, fix_elem_orient);
        data.grid_f = new mfem::GridFunction(data.mesh, mesh_file);
        if (data_type == "fem3d_gf_data_keys")
        {
            mesh_file >> data.keys;
        }
    }
    else if (data_type == "vfem3d_gf_data" || data_type == "vfem3d_gf_data_keys")
    {
        data.field_type = 1;
        data.mesh = new mfem::Mesh(mesh_file, 1, 0, fix_elem_orient);
        data.grid_f = new mfem::GridFunction(data.mesh, mesh_file);
        if (data_type == "vfem3d_gf_data_keys")
        {
            mesh_file >> data.keys;
        }
    }
    else if (data_type == "solution")
    {
        data.mesh = new mfem::Mesh(mesh_file, 1, 0, fix_elem_orient);
        data.grid_f = new mfem::GridFunction(data.mesh, mesh_file);
        data.field_type = (data.grid_f->VectorDim() == 1) ? 0 : 1;
    }
    else if (data_type == "mesh")
    {
        data.mesh = new mfem::Mesh(mesh_file, 1, 0, fix_elem_orient);
        setMeshSolution(data);
        data.field_type = 2;
    }
    else if (data_type == "raw_scalar_2d")
    {
        Array<Array<double> *> vertices;
        Array<Array<int> *> elements;
        Array<int> elem_types;
        string ident;
        int num_patches, num_vert, num_elem, n;
        mesh_file >> std::ws >> ident; // 'patches'
        mesh_file >> num_patches;
        // cout << ident << ' ' << num_patches << endl;
        vertices.SetSize(num_patches);
        vertices = NULL;
        elements.SetSize(num_patches);
        elements = NULL;
        elem_types.SetSize(num_patches);
        elem_types = 0;
        int tot_num_vert = 0;
        int tot_num_elem = 0;
        int mesh_type = 0;
        for (int i = 0; i < num_patches; i++)
        {
            mesh_file >> ws >> ident; // 'vertices'
            mesh_file >> num_vert;
            // cout << '\n' << ident << ' ' << num_vert << endl;
            // read vertices in the format: x y z nx ny nz
            vertices[i] = new Array<double>(6*num_vert);
            Array<double> &verts = *vertices[i];
            for (int j = 0; j < verts.Size(); j++)
            {
                mesh_file >> verts[j];
            }

            mesh_file >> ws >> ident; // 'triangles' or 'quads'
            if (ident == "triangles")
            {
                n = 3, mesh_type |= 1;
            }
            else
            {
                n = 4, mesh_type |= 2;
            }
            elem_types[i] = n;
            mesh_file >> num_elem;
            // cout << ident << ' ' << num_elem << endl;
            elements[i] = new Array<int>(n*num_elem);
            Array<int> &elems = *elements[i];
            for (int j = 0; j < elems.Size(); j++)
            {
                mesh_file >> elems[j];
                elems[j] += tot_num_vert;
            }
            tot_num_vert += num_vert;
            tot_num_elem += num_elem;
        }

        data.mesh = new mfem::Mesh(2, tot_num_vert, tot_num_elem, 0);
        data.sol.SetSize(tot_num_vert);
        data.normals.SetSize(3*tot_num_vert);

        int v_off = 0;
        for (int i = 0; i < num_patches; i++)
        {
            Array<double> &verts = *vertices[i];
            num_vert = verts.Size()/6;
            for (int j = 0; j < num_vert; j++)
            {
                data.mesh->AddVertex(&verts[6*j]);
                data.sol(v_off) = verts[6*j+2];
                data.normals(3*v_off+0) = verts[6*j+3];
                data.normals(3*v_off+1) = verts[6*j+4];
                data.normals(3*v_off+2) = verts[6*j+5];
                v_off++;
            }

            n = elem_types[i];
            Array<int> &elems = *elements[i];
            num_elem = elems.Size()/n;
            // int attr = 1;
            int attr = i + 1;
            if (n == 3)
                for (int j = 0; j < num_elem; j++)
                {
                    data.mesh->AddTriangle(&elems[3*j], attr);
                }
            else
                for (int j = 0; j < num_elem; j++)
                {
                    data.mesh->AddQuad(&elems[4*j], attr);
                }
        }

        if (mesh_type == 1)
        {
            data.mesh->FinalizeTriMesh(1, 0, fix_elem_orient);
        }
        else if (mesh_type == 2)
        {
            data.mesh->FinalizeQuadMesh(1, 0, fix_elem_orient);
        }
        else
        {
            mfem_error("Input data contains mixture of triangles and quads!");
        }

        data.mesh->GenerateBoundaryElements();

        for (int i = num_patches; i > 0; )
        {
            i--;
            delete elements[i];
            delete vertices[i];
        }

        data.field_type = 0;
    }
    else
    {
        data.field_type = -1;
        cerr << "Unknown data format" << endl;
        cerr << data_type << endl;
    }

    if (data.field_type >= 0 && data.field_type <= 2)
    {
        Extrude1DMeshAndSolution(data);
    }
    return data;
}

// Replace a given VectorFiniteElement-based grid function (e.g. from a Nedelec
// or Raviart-Thomas space) with a discontinuous piece-wise polynomial Cartesian
// product vector grid function of the same order.
GridFunction *ProjectVectorFEGridFunction(GridFunction *gf)
{
   if ((gf->VectorDim() == 3) && (gf->FESpace()->GetVDim() == 1))
   {
      int p = gf->FESpace()->GetOrder(0);
      cout << "Switching to order " << p
           << " discontinuous vector grid function..." << endl;
      int dim = gf->FESpace()->GetMesh()->Dimension();
      FiniteElementCollection *d_fec = new L2_FECollection(p, dim, 1);
      FiniteElementSpace *d_fespace =
         new FiniteElementSpace(gf->FESpace()->GetMesh(), d_fec, 3);
      GridFunction *d_gf = new GridFunction(d_fespace);
      d_gf->MakeOwner(d_fec);
      gf->ProjectVectorFieldOn(*d_gf);
      delete gf;
      return d_gf;
   }
   return gf;
}

bool startVisualization(std::string input, std::string data_type, int w, int h) {

    glvis_data = readStream(input, data_type);

    if (glvis_data.field_type < 0 || glvis_data.field_type > 2) {
        return false;
    }

    if (InitVisualization("glvis", 0, 0, w, h)) {
        return false;
    }

    VisualizationSceneScalarData * vs;
    
    double mesh_range = -1.0;
    if (glvis_data.field_type == 0 || glvis_data.field_type == 2)
    {
        if (glvis_data.grid_f)
        {
            glvis_data.grid_f->GetNodalValues(glvis_data.sol);
        }
        if (glvis_data.mesh->SpaceDimension() == 2)
        {
            VisualizationSceneSolution * vss;
            if (glvis_data.field_type == 2)
            {
                paletteSet(4);
            }
            if (glvis_data.normals.Size() > 0)
            {
                vs = vss = new VisualizationSceneSolution(*glvis_data.mesh, glvis_data.sol, &glvis_data.normals);
            }
            else
            {
                vs = vss = new VisualizationSceneSolution(*glvis_data.mesh, glvis_data.sol);
            }
            if (glvis_data.grid_f)
            {
                vss->SetGridFunction(*glvis_data.grid_f);
            }
            if (glvis_data.field_type == 2)
            {
                vs->OrthogonalProjection = 1;
                vs->light = 0;
                vs->Zoom(1.8);
            }
        }
        else if (glvis_data.mesh->SpaceDimension() == 3)
        {
            VisualizationSceneSolution3d * vss;
            vs = vss = new VisualizationSceneSolution3d(*glvis_data.mesh, glvis_data.sol);
            if (glvis_data.grid_f)
            {
                vss->SetGridFunction(glvis_data.grid_f);
            }
            if (glvis_data.field_type == 2)
            {
                if (glvis_data.mesh->Dimension() == 3)
                {
                    paletteSet(4);
                    // paletteSet(11);
                    // Set_Material_And_Light(4,3);
                }
                else
                {
                    paletteSet(4);
                }
                vss->ToggleDrawAxes();
                vss->ToggleDrawMesh();
            }
        }
        if (glvis_data.field_type == 2)
        {
            if (glvis_data.grid_f)
            {
                mesh_range = glvis_data.grid_f->Max() + 1.0;
            }
            else
            {
                mesh_range = glvis_data.sol.Max() + 1.0;
            }
        }
    }
    else if (glvis_data.field_type == 1)
    {
        if (glvis_data.mesh->SpaceDimension() == 2)
        {
            if (glvis_data.grid_f)
            {
                vs = new VisualizationSceneVector(*glvis_data.grid_f);
            }
            else
            {
                vs = new VisualizationSceneVector(*glvis_data.mesh, glvis_data.solu, glvis_data.solv);
            }
        }
        else if (glvis_data.mesh->SpaceDimension() == 3)
        {
            if (glvis_data.grid_f)
            {
                glvis_data.grid_f = ProjectVectorFEGridFunction(glvis_data.grid_f);
                vs = new VisualizationSceneVector3d(*glvis_data.grid_f);
            }
            else
            {
                vs = new VisualizationSceneVector3d(*glvis_data.mesh, glvis_data.solu, glvis_data.solv, glvis_data.solw);
            }
        }
    }

    if (vs)
    {
        // increase the refinement factors if visualizing a GridFunction
        if (glvis_data.grid_f)
        {
            vs->AutoRefine();
            vs->SetShading(2, true);
        }
        if (mesh_range > 0.0)
        {
            vs->SetValueRange(-mesh_range, mesh_range);
            vs->SetAutoscale(0);
        }
        if (glvis_data.mesh->SpaceDimension() == 2 && glvis_data.field_type == 2)
        {
            SetVisualizationScene(vs, 2, glvis_data.keys.c_str());
        }
        else
        {
            SetVisualizationScene(vs, 3, glvis_data.keys.c_str());
        }
    }
    SendExposeEvent();
    return true;
}

void iterVisualization() {
    GetAppWindow()->mainIter();
}

void disableKeyHandling() {
    SDL_EventState(SDL_KEYDOWN, SDL_DISABLE);
    SDL_EventState(SDL_KEYUP, SDL_DISABLE);
}

void enableKeyHandling() {
    SDL_EventState(SDL_KEYDOWN, SDL_ENABLE);
    SDL_EventState(SDL_KEYUP, SDL_ENABLE);
}

void setKeyboardListeningElementId(const std::string & elem_id) {
  SDL_SetHint(SDL_HINT_EMSCRIPTEN_KEYBOARD_ELEMENT, elem_id.c_str());
}

class SceneModel {
    VisualizationSceneScalarData * _scene;
    SceneModel(VisualizationSceneScalarData * scene)
        : _scene(scene) { }

public:
    static SceneModel create() {
        return SceneModel(dynamic_cast<VisualizationSceneScalarData*>
                                        (GetVisualizationScene()));
    }

    int getDrawAxes() const { return _scene->GetDrawAxes(); }
    int getShading() const { return _scene->GetShading(); }
    int getColorbar() const { return _scene->GetColorbar(); }
    int getRuler() const { return _scene->GetRuler(); }

    void setDrawAxes (int axes) {
        _scene->SetDrawAxes(axes);
        SendExposeEvent();
    }
    void setShading (int shade) {
        _scene->SetShading(shade, true);
        SendExposeEvent();
    }
    void setColorbar (int colorbar) { _scene->SetColorbar(colorbar); }
    void setRuler (int ruler) { _scene->SetRuler(ruler); }

};

} // namespace js

namespace em = emscripten;
EMSCRIPTEN_BINDINGS(js_funcs) {
    em::function("startVisualization", &js::startVisualization);
    em::function("iterVisualization", &js::iterVisualization);
    em::function("sendExposeEvent", &SendExposeEvent);
    em::function("disableKeyHanding", &js::disableKeyHandling);
    em::function("enableKeyHandling", &js::enableKeyHandling);
    em::function("setKeyboardListeningElementId", js::setKeyboardListeningElementId);
    em::function("getTextureMode", &GetUseTexture);
    em::function("setTextureMode", &SetUseTexture);
    em::function("resizeWindow", &ResizeWindow);
}
EMSCRIPTEN_BINDINGS(model) {
    em::class_<js::SceneModel>("SceneModel")
        .property("shading", &js::SceneModel::getShading, &js::SceneModel::setShading)
        .property("colorbar", &js::SceneModel::getColorbar, &js::SceneModel::setColorbar)
        .property("drawaxes", &js::SceneModel::getDrawAxes, &js::SceneModel::setDrawAxes)
        .property("ruler", &js::SceneModel::getRuler, &js::SceneModel::setRuler)
        .class_function("create", &js::SceneModel::create);
}
