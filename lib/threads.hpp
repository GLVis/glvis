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

#ifndef GLVIS_THREADS
#define GLVIS_THREADS

#include <pthread.h>

class GLVisCommand
{
private:
   // Pointers to global GLVis data
   VisualizationSceneScalarData **vs;
   Mesh          **mesh;
   GridFunction  **grid_f;
   Vector         *sol;
   int            *autoscale, *keep_attr;

   pthread_mutex_t glvis_mutex;
   pthread_cond_t  glvis_cond;
   int num_waiting;
   bool terminating;
   int pfd[2];  // pfd[0] -- reading, pfd[1] -- writing

   enum { NO_COMMAND = 0, NEW_MESH_AND_SOLUTION = 1 };

   // command to be executed
   int command;

   // command arguments
   Mesh         *new_m;
   GridFunction *new_g;

   int lock();
   int signal();
   void unlock();

public:
   // called by the main execution thread
   GLVisCommand(VisualizationSceneScalarData **_vs, Mesh **_mesh,
                GridFunction **_grid_f, Vector *_sol, int *_autoscale,
                int *_keep_attr);

   // to be used by the main execution (visualization) thread
   int ReadFD() { return pfd[0]; }

   // to be used worker threads
   int KeepAttrib() { return *keep_attr; } // may need to sync this

   // called by worker threads
   int NewMeshAndSolution(Mesh *_new_m, GridFunction *_new_g);

   // called by the main execution thread
   int Execute();

   // called by the main execution thread
   void Terminate();

   // called by the main execution thread
   ~GLVisCommand();
};

extern GLVisCommand *glvis_command;

class communication_thread
{
private:
   // streams to read data from
   Array<istream *> &is;

   // data that may be dynamically allocated by the thread
   Mesh *new_m;
   GridFunction *new_g;
   string ident;

   // thread id
   pthread_t tid;

   static void cancel_off()
   {
      pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, NULL);
   }

   static void cancel_on()
   {
      pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
   }

   static void *execute(void *);

public:
   communication_thread(Array<istream *> &_is);

   ~communication_thread();
};

#endif
