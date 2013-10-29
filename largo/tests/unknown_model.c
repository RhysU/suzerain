//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// largo 0.0.1: largo - slow growth terms for turbulence simulations
// http://pecos.ices.utexas.edu/
//
// Copyright (C) 2011, 2012, 2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
// $Id$

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <largo/largo.h>

// Include FCTX and silence useless warnings
#ifdef __INTEL_COMPILER
#pragma warning(push,disable:981)
#endif
#include "fct.h"
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif

FCT_BGN()
{
    FCT_QTEST_BGN(something)
    {
        // Produce non-NULL address
        char dummy = 'c';
        void * generic_workspace = &dummy;
        fct_chk(generic_workspace != NULL);

        // Allocate workspace with unknown model name produces NULL pointer
        largo_allocate (&generic_workspace, "unknown_model", 7, 2, 0, "dns");
        fct_chk(generic_workspace == NULL);

        // Deallocation of the NULL pointer is okay
        largo_deallocate (&generic_workspace);
        fct_chk(generic_workspace == NULL);
    }
    FCT_QTEST_END();
}
FCT_END()
