#ifndef __CGLCONICCONFIG_H__
#define __CGLCONICCONFIG_H__

#ifdef HAVE_CONFIG_H
#ifdef CGLCONIC_BUILD
#include "config.h"
#else
#include "config_cglconic.h"
#endif

#else /* HAVE_CONFIG_H */

#ifdef CGLCONIC_BUILD
#include "config_default.h"
#else
#include "config_cglconic_default.h"
#endif

#endif /* HAVE_CONFIG_H */

#endif /*__CGLCONICCONFIG_H__*/
