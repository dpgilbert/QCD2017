#ifndef QCDCHECK_H
#define QCDCHECK_H

#include <iostream>
#include <vector>
#include <map>
#include <cmath>

#include "TH1.h"

#include "../CORE/Tools/utils.h"
#include "../CORE/Tools/goodrun.h"
#include "../CORE/Tools/dorky/dorky.h"
#include "../CORE/Tools/badEventFilter.h"

#include "../MT2CORE/mt2tree.h"
#include "../MT2CORE/sigSelections.h"
#include "../MT2CORE/SR.h"

class QCDcheck {
 public:
  QCDcheck();
  ~QCDcheck();

  void loop (const TChain* chain, const char* sample, const char* outdir, const int maxEvt);
};

#endif
