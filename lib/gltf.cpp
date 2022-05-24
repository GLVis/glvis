// Copyright (c) 2010-2022, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-443271.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see https://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "gltf.hpp"
#include "aux_vis.hpp" // SaveAsPNG

using namespace std;


const char *glTF_Builder::tensorTypes[] =
{
   "SCALAR", "VEC2", "VEC3", "VEC4", "MAT2", "MAT3", "MAT4"
};

glTF_Builder::buffer_id
glTF_Builder::addBuffer(const string &bufferName)
{
   buffers.resize(buffers.size() + 1);
   auto &buf = buffers.back();
   buf.uri.value = file_prefix + "." + bufferName + ".bin";
   buf.uri.valid = true;
   buf.byteLength.value = 0;
   buf.byteLength.valid = true;
   buf.file.reset(new ofstream(buf.uri.value, ios::out | ios::binary));

   return {(unsigned)buffers.size() - 1};
}

glTF_Builder::buffer_view_id
glTF_Builder::addBufferView(buffer_id buffer,
                            const void *data,
                            size_t byteLength,
                            size_t byteStride,
                            size_t byteAlign,
                            target_type target)
{
   if (buffer.id >= buffers.size()) { return {INVALID_ID}; }

   buffer_views.resize(buffer_views.size() + 1);
   auto &buf_view = buffer_views.back();
   auto &buf = buffers[buffer.id];

   buf_view.buffer.value = buffer.id;
   buf_view.buffer.valid = true;

   const unsigned buf_offset = buf.byteLength.value;
   const unsigned new_offset = byteAlign*((buf_offset+byteAlign-1)/byteAlign);
   buf_view.byteOffset.value = new_offset;
   buf_view.byteOffset.valid = true;

   buf_view.byteLength.value = byteLength;
   buf_view.byteLength.valid = true;

   if (target == target_type::ARRAY_BUFFER)
   {
      buf_view.byteStride.value = byteStride;
      buf_view.byteStride.valid = true;
   }

   buf_view.target.value = (unsigned)target;
   buf_view.target.valid = true;

   // append padding to file
   for (unsigned i = buf_offset; i != new_offset; ++i) { buf.file->put('\0'); }
   // write data to file
   buf.file->write(reinterpret_cast<const char *>(data), byteLength);

   buf.byteLength.value = new_offset + byteLength;

   return {(unsigned)buffer_views.size() - 1};
}

void glTF_Builder::appendToBufferView(buffer_view_id bufferView,
                                      const void *data,
                                      size_t byteLength)
{
   if (bufferView.id >= buffer_views.size()) { return; }

   auto &buf_view = buffer_views[bufferView.id];
   auto &buf = buffers[buf_view.buffer.value];

   buf_view.byteLength.value += byteLength;

   buf.file->write(reinterpret_cast<const char *>(data), byteLength);
   buf.byteLength.value += byteLength;
}

glTF_Builder::accessor_id
glTF_Builder::addAccessor(buffer_view_id bufferView,
                          size_t byteOffset,
                          component_type componentType,
                          size_t count,
                          tensor_type tensorType)
{
   if (bufferView.id >= buffer_views.size() || count == 0)
   {
      return {INVALID_ID};
   }

   accessors.resize(accessors.size() + 1);
   auto &acc = accessors.back();

   acc.bufferView.value = bufferView.id;
   acc.bufferView.valid = true;

   acc.byteOffset.value = byteOffset;
   acc.byteOffset.valid = true;

   acc.componentType.value = (unsigned)componentType;
   acc.componentType.valid = true;

   acc.count.value = count;
   acc.count.valid = true;

   acc.type.value = tensorTypes[(unsigned)tensorType];
   acc.type.valid = true;

   // Note: acc.min and acc.max remain invalid and will not be written too file

   if (componentType != component_type::FLOAT &&
       buffer_views[bufferView.id].target.value !=
       (unsigned)target_type::ELEMENT_ARRAY_BUFFER)
   {
      acc.normalized.value = true;
      acc.normalized.valid = true;
   }

   return {(unsigned)accessors.size() - 1};
}

glTF_Builder::accessor_id
glTF_Builder::addAccessorVec2f(buffer_view_id bufferView,
                               size_t byteOffset,
                               size_t count,
                               vec2f min,
                               vec2f max)
{
   auto id = addAccessor(bufferView,
                         byteOffset,
                         component_type::FLOAT,
                         count,
                         tensor_type::VEC2);

   if (id.id != INVALID_ID)
   {
      auto &acc = accessors[id.id];
      acc.min.value.assign(min.begin(), min.end());
      acc.min.valid = true;
      acc.max.value.assign(max.begin(), max.end());
      acc.max.valid = true;
   }

   return id;
}

glTF_Builder::accessor_id
glTF_Builder::addAccessorVec3f(buffer_view_id bufferView,
                               size_t byteOffset,
                               size_t count,
                               vec3f min,
                               vec3f max)
{
   auto id = addAccessor(bufferView,
                         byteOffset,
                         component_type::FLOAT,
                         count,
                         tensor_type::VEC3);

   if (id.id != INVALID_ID)
   {
      auto &acc = accessors[id.id];
      acc.min.value.assign(min.begin(), min.end());
      acc.min.valid = true;
      acc.max.value.assign(max.begin(), max.end());
      acc.max.valid = true;
   }

   return id;
}

glTF_Builder::image_id
glTF_Builder::addImage(const string &imageName,
                       int width,
                       int height,
                       const color4f *pixels)
{
#ifndef GLVIS_USE_LIBPNG

   return {INVALID_ID};

#else

   images.resize(images.size() + 1);
   auto &img = images.back();

   img.uri.value = file_prefix + "." + imageName + ".png";
   img.uri.valid = true;

   img.name.value = imageName;
   img.name.valid = true;

   // write the image
   auto get_row = [&](int row, void *pxls)
   {
      auto pxls_out = reinterpret_cast<array<uint8_t,4>*>(pxls);
      auto pxls_in = pixels + row*width;
      for (int i = 0; i < width; ++i)
      {
         for (int j = 0; j < 4; ++j)
         {
            pxls_out[i][j] = std::min(int(pxls_in[i][j]*256), 255);
         }
      }
   };
   SaveAsPNG(img.uri.value.c_str(), width, height,
             /* is_hidpi: */ false, /* with_alpha: */ true, get_row);

   return {(unsigned)images.size() - 1};

#endif // GLVIS_USE_LIBPNG
}

glTF_Builder::sampler_id
glTF_Builder::addSampler(mag_filter magFilter,
                         min_filter minFilter,
                         wrap_type wrapS,
                         wrap_type wrapT)
{
   samplers.resize(samplers.size() + 1);
   auto &sampler = samplers.back();

   sampler.magFilter.value = (unsigned)magFilter;
   sampler.magFilter.valid = true;

   sampler.minFilter.value = (unsigned)minFilter;
   sampler.minFilter.valid= true;

   sampler.wrapS.value = (unsigned)wrapS;
   sampler.wrapS.valid = true;

   sampler.wrapT.value = (unsigned)wrapT;
   sampler.wrapT.valid = true;

   return {(unsigned)samplers.size() - 1};
}

glTF_Builder::texture_id
glTF_Builder::addTexture(sampler_id sampler, image_id source)
{
   if (sampler.id >= samplers.size() || source.id >= images.size())
   {
      return {INVALID_ID};
   }

   textures.resize(textures.size() + 1);
   auto &tex = textures.back();

   tex.sampler.value = sampler.id;
   tex.sampler.valid = true;

   tex.source.value = source.id;
   tex.source.valid = true;

   return {(unsigned)textures.size() - 1};
}

glTF_Builder::material_id
glTF_Builder::addMaterial(const string &materialName,
                          const pbr_matallic_roughness &pbrMetallicRoughness,
                          bool doubleSided)
{
   if (pbrMetallicRoughness.haveTexture &&
       pbrMetallicRoughness.baseColorTexture.id >= textures.size())
   {
      return {INVALID_ID};
   }

   materials.resize(materials.size() + 1);
   auto &mat = materials.back();

   mat.name.value = materialName;
   mat.name.valid = true;

   auto &pbr = mat.pbrMetallicRoughness.value;

   pbr.baseColorFactor.value = pbrMetallicRoughness.baseColorFactor;
   pbr.baseColorFactor.valid = true;

   auto &tex_info = pbr.baseColorTexture.value;

   tex_info.index.value = pbrMetallicRoughness.baseColorTexture.id;
   tex_info.index.valid = true;

   tex_info.texCoord.value = 0;
   tex_info.texCoord.valid = true;

   pbr.baseColorTexture.valid = pbrMetallicRoughness.haveTexture;

   pbr.metallicFactor.value = pbrMetallicRoughness.metallicFactor;
   pbr.metallicFactor.valid = true;

   pbr.roughnessFactor.value = pbrMetallicRoughness.roughnessFactor;
   pbr.roughnessFactor.valid = true;

   mat.pbrMetallicRoughness.valid = true;

   mat.doubleSided.value = doubleSided;
   mat.doubleSided.valid = true;

   return {(unsigned)materials.size() - 1};
}

glTF_Builder::mesh_id
glTF_Builder::addMesh(const string &meshName)
{
   meshes.resize(meshes.size() + 1);
   auto &mesh = meshes.back();

   mesh.name.value = meshName;
   mesh.name.valid = true;

   return {(unsigned)meshes.size() - 1};
}

void glTF_Builder::addMeshTriangles(mesh_id mesh,
                                    accessor_id vertexPositions,
                                    accessor_id vertexNormals,
                                    accessor_id vertexTexCoords0,
                                    accessor_id vertexIndices,
                                    material_id material)
{
   if (mesh.id >= meshes.size() || vertexPositions.id >= accessors.size())
   { return; }

   auto &primitives = meshes[mesh.id].primitives;
   primitives.resize(primitives.size() + 1);
   auto &pri = primitives.back();

   pri.attributes.value.POSITION.value = vertexPositions.id;
   pri.attributes.value.POSITION.valid = true;

   if (vertexNormals.id < accessors.size())
   {
      pri.attributes.value.NORMAL.value = vertexNormals.id;
      pri.attributes.value.NORMAL.valid = true;
   }

   if (vertexTexCoords0.id < accessors.size())
   {
      pri.attributes.value.TEXCOORD_0.value = vertexTexCoords0.id;
      pri.attributes.value.TEXCOORD_0.valid = true;
   }

   pri.attributes.valid = true;

   if (vertexIndices.id < accessors.size())
   {
      pri.indices.value = vertexIndices.id;
      pri.indices.valid = true;
   }

   if (material.id < materials.size())
   {
      pri.material.value = material.id;
      pri.material.valid = true;
   }

   // pri.mode remains undefined since default is 4 = TRIANGLES
}

void glTF_Builder::addMeshLines(mesh_id mesh,
                                accessor_id vertexPositions,
                                accessor_id vertexTexcoords0,
                                accessor_id vertexColors0,
                                material_id material)
{
   if (mesh.id >= meshes.size()) { return; }

   auto &primitives = meshes[mesh.id].primitives;
   primitives.resize(primitives.size() + 1);
   auto &pri = primitives.back();

   pri.attributes.value.POSITION.value = vertexPositions.id;
   pri.attributes.value.POSITION.valid = true;

   if (vertexTexcoords0.id < accessors.size())
   {
      pri.attributes.value.TEXCOORD_0.value = vertexTexcoords0.id;
      pri.attributes.value.TEXCOORD_0.valid = true;
   }
   else if (vertexColors0.id < accessors.size())
   {
      pri.attributes.value.COLOR_0.value = vertexColors0.id;
      pri.attributes.value.COLOR_0.valid = true;
   }

   // NORMAL remains undefined

   pri.attributes.valid = true;

   // pri.indices remain undefined

   if (material.id < materials.size())
   {
      pri.material.value = material.id;
      pri.material.valid = true;
   }

   pri.mode.value = 1; // = LINES
   pri.mode.valid = true;
}

glTF_Builder::node_id
glTF_Builder::addNode(const string &nodeName)
{
   nodes.resize(nodes.size() + 1);
   auto &node = nodes.back();

   node.name.value = nodeName;
   node.name.valid = true;

   return {(unsigned)nodes.size() - 1};
}

void glTF_Builder::addNodeMesh(node_id node, mesh_id mesh)
{
   if (node.id >= nodes.size()) { return; }

   nodes[node.id].mesh.value = mesh.id;
   nodes[node.id].mesh.valid = true;
}

void glTF_Builder::addNodeScale(node_id node, vec3f scale)
{
   if (node.id >= nodes.size()) { return; }

   nodes[node.id].scale.value = scale;
   nodes[node.id].scale.valid = true;
}

void glTF_Builder::addNodeTranslation(node_id node, vec3f translation)
{
   if (node.id >= nodes.size()) { return; }

   nodes[node.id].translation.value = translation;
   nodes[node.id].translation.valid = true;
}

void glTF_Builder::getMaterialPBRMR(material_id material,
                                    pbr_matallic_roughness &pbr_mr_copy)
{
   if (material.id >= materials.size()) { return; }

   auto &mat = materials[material.id];
   auto &pbr = mat.pbrMetallicRoughness.value;

   pbr_mr_copy.haveTexture = pbr.baseColorTexture.valid;
   pbr_mr_copy.baseColorFactor = pbr.baseColorFactor.value;
   pbr_mr_copy.baseColorTexture = {pbr.baseColorTexture.value.index.value};
   pbr_mr_copy.metallicFactor = pbr.metallicFactor.value;
   pbr_mr_copy.roughnessFactor = pbr.roughnessFactor.value;
}

int glTF_Builder::writeFile()
{
   if (nodes.size() == 0)
   {
      return 1;
   }

   ofstream gltf(file_prefix + ".gltf");
   gltf.precision(8);
   gltf.setf(ios::boolalpha);

   // ~~~ Tutorial ~~~
   // https://github.com/KhronosGroup/glTF-Tutorials/blob/master/gltfTutorial/README.md

   // https://www.khronos.org/registry/glTF/specs/2.0/glTF-2.0.html#reference-scene
   gltf <<
        "{\n"
        "  \"scene\": 0,\n"
        "  \"scenes\" : [ {\n"
        "      \"nodes\" : [";
   for (size_t i = 0; i != nodes.size(); ++i) { gltf << sep(i) << ' ' << i; }
   gltf <<
        " ]\n"
        "  } ],\n\n";

   // https://www.khronos.org/registry/glTF/specs/2.0/glTF-2.0.html#reference-node
   gltf << "  \"nodes\" : [";
   for (size_t i = 0; i != nodes.size(); ++i)
   {
      gltf << sep(i) << " {";
      int pos = 0;
      print_node(gltf, pos, "\n    ", nodes[i].name);
      print_node(gltf, pos, "\n    ", nodes[i].mesh);
      print_node(gltf, pos, "\n    ", nodes[i].scale);
      print_node(gltf, pos, "\n    ", nodes[i].translation);
      gltf << "\n  }";
   }
   gltf << " ],\n\n";

   // https://www.khronos.org/registry/glTF/specs/2.0/glTF-2.0.html#reference-mesh
   gltf << "  \"meshes\" : [";
   for (size_t i = 0; i != meshes.size(); ++i)
   {
      gltf << sep(i) << " {";
      int pos = 0;
      print_node(gltf, pos, "\n    ", meshes[i].name);
      gltf << sep(pos++) << "\n    \"primitives\" : [";
      auto &primitives = meshes[i].primitives;
      for (size_t j = 0; j != primitives.size(); ++j)
      {
         gltf << sep(j) << " {";
         int pos2 = 0;
         print_node(gltf, pos2, "\n      ", primitives[j].attributes);
         print_node(gltf, pos2, "\n      ", primitives[j].indices);
         print_node(gltf, pos2, "\n      ", primitives[j].material);
         print_node(gltf, pos2, "\n      ", primitives[j].mode);
         gltf << "\n    }";
      }
      gltf << " ]\n  }";
   }
   gltf << " ],\n\n";

   // https://www.khronos.org/registry/glTF/specs/2.0/glTF-2.0.html#reference-material
   gltf << "  \"materials\" : [";
   for (size_t i = 0; i != materials.size(); ++i)
   {
      gltf << sep(i) << " {";
      int pos = 0;
      print_node(gltf, pos, "\n    ", materials[i].name);
      print_node(gltf, pos, "\n    ", materials[i].pbrMetallicRoughness);
      print_node(gltf, pos, "\n    ", materials[i].doubleSided);
      gltf << "\n  }";
   }
   gltf << " ],\n\n";

   gltf << "  \"textures\" : [";
   for (size_t i = 0; i != textures.size(); ++i)
   {
      gltf << sep(i) << " {";
      int pos = 0;
      print_node(gltf, pos, "\n    ", textures[i].sampler);
      print_node(gltf, pos, "\n    ", textures[i].source);
      gltf << "\n  }";
   }
   gltf << " ],\n\n";

   gltf << "  \"images\" : [";
   for (size_t i = 0; i != images.size(); ++i)
   {
      gltf << sep(i) << " {";
      int pos = 0;
      print_node(gltf, pos, "\n    ", images[i].name);
      print_node(gltf, pos, "\n    ", images[i].uri);
      gltf << "\n  }";
   }
   gltf << " ],\n\n";

   // https://www.khronos.org/registry/glTF/specs/2.0/glTF-2.0.html#reference-sampler
   // see also: PaletteState::{ToTextureDiscrete(),ToTextureSmooth()}
   gltf << "  \"samplers\" : [";
   for (size_t i = 0; i != samplers.size(); ++i)
   {
      gltf << sep(i) << " {";
      int pos = 0;
      print_node(gltf, pos, "\n    ", samplers[i].magFilter);
      print_node(gltf, pos, "\n    ", samplers[i].minFilter);
      print_node(gltf, pos, "\n    ", samplers[i].wrapS);
      print_node(gltf, pos, "\n    ", samplers[i].wrapT);
      gltf << "\n  }";
   }
   gltf << " ],\n\n";

   // https://www.khronos.org/registry/glTF/specs/2.0/glTF-2.0.html#reference-buffer
   gltf << "  \"buffers\" : [";
   for (size_t i = 0; i != buffers.size(); ++i)
   {
      gltf << sep(i) << " {";
      int pos = 0;
      print_node(gltf, pos, "\n    ", buffers[i].uri);
      print_node(gltf, pos, "\n    ", buffers[i].byteLength);
      gltf << "\n  }";
   }
   gltf << " ],\n\n";

   // https://www.khronos.org/registry/glTF/specs/2.0/glTF-2.0.html#reference-bufferview
   gltf << "  \"bufferViews\" : [";
   for (size_t i = 0; i != buffer_views.size(); ++i)
   {
      gltf << sep(i) << " {";
      int pos = 0;
      print_node(gltf, pos, "\n    ", buffer_views[i].buffer);
      print_node(gltf, pos, "\n    ", buffer_views[i].byteOffset);
      print_node(gltf, pos, "\n    ", buffer_views[i].byteLength);
      print_node(gltf, pos, "\n    ", buffer_views[i].byteStride);
      print_node(gltf, pos, "\n    ", buffer_views[i].target);
      gltf << "\n  }";
   }
   gltf << " ],\n\n";

   // https://www.khronos.org/registry/glTF/specs/2.0/glTF-2.0.html#reference-accessor
   gltf << "  \"accessors\" : [";
   for (size_t i = 0; i != accessors.size(); ++i)
   {
      gltf << sep(i) << " {";
      int pos = 0;
      print_node(gltf, pos, "\n    ", accessors[i].bufferView);
      print_node(gltf, pos, "\n    ", accessors[i].byteOffset);
      print_node(gltf, pos, "\n    ", accessors[i].componentType);
      print_node(gltf, pos, "\n    ", accessors[i].count);
      print_node(gltf, pos, "\n    ", accessors[i].type);
      print_node(gltf, pos, "\n    ", accessors[i].min);
      print_node(gltf, pos, "\n    ", accessors[i].max);
      print_node(gltf, pos, "\n    ", accessors[i].normalized);
      gltf << "\n  }";
   }
   gltf << " ],\n\n";

   // https://www.khronos.org/registry/glTF/specs/2.0/glTF-2.0.html#reference-asset
   gltf <<
        "  \"asset\" : {\n"
        "    \"version\" : \"2.0\",\n"
        "    \"generator\" : \"GLVis\"\n"
        "  }\n"
        "}\n";

   return 0;
}
