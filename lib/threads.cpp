// Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443271. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see http://glvis.googlecode.com.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.

#include <unistd.h>    // pipe, fcntl, write
#include <fcntl.h>     // fcntl
#include <errno.h>     // errno, EAGAIN
#include <stdio.h>     // perror
#include "visual.hpp"

GLVisCommand *glvis_command = NULL;

GLVisCommand::GLVisCommand(
   VisualizationSceneScalarData **_vs, Mesh **_mesh, GridFunction **_grid_f,
   Vector *_sol, int *_autoscale, int *_keep_attr)
{
   vs        = _vs;
   mesh      = _mesh;
   grid_f    = _grid_f;
   sol       = _sol;
   autoscale = _autoscale;
   keep_attr = _keep_attr;

   pthread_mutex_init(&glvis_mutex, NULL);
   pthread_cond_init(&glvis_cond, NULL);
   num_waiting = 0;
   terminating = false;
   if (pipe(pfd) == -1)
   {
      perror("pipe()");
      exit(EXIT_FAILURE);
   }
   int flag = fcntl(pfd[0], F_GETFL);
   fcntl(pfd[0], F_SETFL, flag | O_NONBLOCK);

   command = NO_COMMAND;
}

int GLVisCommand::lock()
{
   int my_id;
   pthread_mutex_lock(&glvis_mutex);
   if (terminating)
   {
      pthread_mutex_unlock(&glvis_mutex);
      return -1;
   }
   my_id = num_waiting++;
   while (my_id > 0)
   {
      pthread_cond_wait(&glvis_cond, &glvis_mutex);
      if (terminating)
      {
         num_waiting--;
         pthread_mutex_unlock(&glvis_mutex);
         return -1;
      }
      my_id--;
   }
   pthread_mutex_unlock(&glvis_mutex);
   return 0;
}

int GLVisCommand::signal()
{
   char c = 's';
   if (write(pfd[1], &c, 1) != 1)
      return -1;
   return 0;
}

void GLVisCommand::unlock()
{
   pthread_mutex_lock(&glvis_mutex);
   num_waiting--;
   if (num_waiting > 0)
      pthread_cond_broadcast(&glvis_cond);
   pthread_mutex_unlock(&glvis_mutex);
}

int GLVisCommand::NewMeshAndSolution(Mesh *_new_m, GridFunction *_new_g)
{
   if (lock() < 0)
      return -1;
   command = NEW_MESH_AND_SOLUTION;
   new_m = _new_m;
   new_g = _new_g;
   if (signal() < 0)
      return -2;
   return 0;
}

extern GridFunction *ProjectVectorFEGridFunction(GridFunction*);

int GLVisCommand::Execute()
{
   char c;
   int n = read(pfd[0], &c, 1);

   if (n == -1 && errno == EAGAIN)
      return 1;
   if (n != 1 || c != 's')
      return -1;

   switch (command)
   {
   case NO_COMMAND:
      break;

   case NEW_MESH_AND_SOLUTION:
      if (new_m->Dimension() == (*mesh)->Dimension() &&
          new_g->VectorDim() == (*grid_f)->VectorDim())
      {
         if (new_m->Dimension() == 2)
         {
            if (new_g->VectorDim() == 1)
            {
               VisualizationSceneSolution *vss =
                  dynamic_cast<VisualizationSceneSolution *>(*vs);
               new_g->GetNodalValues(*sol);
               vss->NewMeshAndSolution(new_m, sol, new_g, *autoscale);
            }
            else
            {
               VisualizationSceneVector *vsv =
                  dynamic_cast<VisualizationSceneVector *>(*vs);
               vsv->NewMeshAndSolution(*new_g, *autoscale);
            }
         }
         else
         {
            if (new_g->VectorDim() == 1)
            {
               VisualizationSceneSolution3d *vss =
                  dynamic_cast<VisualizationSceneSolution3d *>(*vs);
               new_g->GetNodalValues(*sol);
               vss->NewMeshAndSolution(new_m, sol, new_g, *autoscale);
            }
            else
            {
               new_g = ProjectVectorFEGridFunction(new_g);
               VisualizationSceneVector3d *vss =
                  dynamic_cast<VisualizationSceneVector3d *>(*vs);
               vss->NewMeshAndSolution(new_m, new_g, *autoscale);
            }
         }
         delete (*grid_f);
         *grid_f = new_g;
         delete (*mesh);
         *mesh = new_m;

         (*vs)->Draw();
      }
      else
      {
         cout << "Stream: field type does not match!" << endl;
         delete new_g;
         delete new_m;
      }
      break;
   }

   command = NO_COMMAND;
   unlock();
   return 0;
}

void GLVisCommand::Terminate()
{
   char c;
   int n = read(pfd[0], &c, 1);

   pthread_mutex_lock(&glvis_mutex);
   terminating = true;
   pthread_mutex_unlock(&glvis_mutex);
   if (n == 1 && c == 's')
   {
      switch (command)
      {
      case NEW_MESH_AND_SOLUTION:
         delete new_g;
         delete new_m;
         break;
      }
      unlock();
   }
   pthread_mutex_lock(&glvis_mutex);
   if (num_waiting > 0)
      pthread_cond_broadcast(&glvis_cond);
   pthread_mutex_unlock(&glvis_mutex);
}

GLVisCommand::~GLVisCommand()
{
   if (num_waiting > 0)
      cout << "\nGLVisCommand::~GLVisCommand() : num_waiting = "
           << num_waiting << '\n' << endl;
   close(pfd[0]);
   close(pfd[1]);
   pthread_cond_destroy(&glvis_cond);
   pthread_mutex_destroy(&glvis_mutex);
}

communication_thread::communication_thread(Array<istream *> &_is)
   : is(_is)
{
   new_m = NULL;
   new_g = NULL;

   if (is.Size() > 0)
      pthread_create(&tid, NULL, communication_thread::execute, this);
}

communication_thread::~communication_thread()
{
   if (is.Size() > 0)
   {
      pthread_cancel(tid);
      pthread_join(tid, NULL);
   }

   delete new_g;
   delete new_m;
}

void *communication_thread::execute(void *p)
{
   communication_thread *_this = (communication_thread *)p;

   while (1)
   {
      *_this->is[0] >> ws;

      _this->cancel_off();

      *_this->is[0] >> _this->ident;
      if (!(*_this->is[0]))
         break;

      if (_this->ident == "solution" || _this->ident == "parallel")
      {
         if (_this->ident == "solution")
         {
            _this->new_m = new Mesh(*_this->is[0], 1, 0);
            if (!(*_this->is[0]))
               break;
            _this->new_g = new GridFunction(_this->new_m, *_this->is[0]);
            if (!(*_this->is[0]))
               break;
         }
         else if (_this->ident == "parallel")
         {
            Array<Mesh *> mesh_array;
            Array<GridFunction *> gf_array;
            int proc, nproc, np = 0;
            int keep_attr = glvis_command->KeepAttrib();
            do
            {
               istream &isock = *_this->is[np];
               isock >> nproc >> proc >> ws;
               isock >> _this->ident >> ws; // "solution"
               mesh_array.SetSize(nproc);
               gf_array.SetSize(nproc);
               mesh_array[proc] = new Mesh(isock, 1, 0);
               if (!keep_attr)
               {
                  // set element and boundary attributes to proc+1
                  for (int i = 0; i < mesh_array[proc]->GetNE(); i++)
                     mesh_array[proc]->GetElement(i)->SetAttribute(proc+1);
                  for (int i = 0; i < mesh_array[proc]->GetNBE(); i++)
                     mesh_array[proc]->GetBdrElement(i)->SetAttribute(proc+1);
               }
               gf_array[proc] = new GridFunction(mesh_array[proc], isock);
               np++;
               if (np == nproc)
                  break;
               *_this->is[np] >> _this->ident >> ws; // "parallel"
            }
            while (1);
            _this->new_m = new Mesh(mesh_array, nproc);
            _this->new_g = new GridFunction(_this->new_m, gf_array, nproc);

            for (int p = 0; p < nproc; p++)
            {
               delete gf_array[nproc-1-p];
               delete mesh_array[nproc-1-p];
            }
            gf_array.DeleteAll();
            mesh_array.DeleteAll();
         }

         // cout << "Stream: new solution" << endl;

         if (glvis_command->NewMeshAndSolution(_this->new_m, _this->new_g))
            goto comm_terminate;

         _this->new_m = NULL;
         _this->new_g = NULL;
      }
      else
      {
         cout << "Stream: unknown command: " << _this->ident << endl;
      }

      _this->cancel_on();
   }

   cout << "Stream: end of input." << endl;

comm_terminate:
   for (int i = 0; i < _this->is.Size(); i++)
   {
      socketstream *isock = dynamic_cast<socketstream *>(_this->is[i]);
      if (isock)
         isock->close();
   }
   _this->cancel_on();

   return p;
}
