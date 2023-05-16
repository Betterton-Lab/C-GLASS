#include "cglass/cell_list.hpp"

double CellList::_min_cell_length_ = -1;
int CellList::_n_dim_ = -1;
int CellList::_n_periodic_ = -1;
int CellList::_n_cells_1d_ = -1;
double CellList::_cell_length_ = -1;
bool CellList::_no_init_ = true;

void CellList::SetMinCellLength(double l) {
  if (l > _min_cell_length_) {
    Logger::Debug("Setting minimum cell length to be %2.2f", l);
    //_min_cell_length_ = 3.75580*1.5; //l;
    _min_cell_length_ = l;
  }
  
}
double CellList::GetCellLength() { return _cell_length_; }

void CellList::Init(int n_dim, int n_periodic, double system_radius) {
  _n_cells_1d_ = 1; //(int)floor(2 * system_radius / _min_cell_length_);
  //_n_cells_1d_ = (int)floor(240 / _min_cell_length_);
  Logger::Info("Number of cells %i", _n_cells_1d_);
  _no_init_ = false;
#ifdef TRACE
  if (_n_cells_1d_ > 20) {
    Logger::Warning("Simulation run in trace mode with a large number of "
                    "cells in cell list (%d).",
                    _n_cells_1d_);
    fprintf(stderr, "Continue anyway? (y/N) ");
    char c;
    if (std::cin.peek() != 'y') {
      Logger::Error("Terminating simulation by user request");
    } else if (!(std::cin >> c)) {
      Logger::Error("Invalid input");
    } else {
      fprintf(stderr, "Resuming simulation\n");
      std::cin.ignore();
    }
  }
#endif
  if (_n_cells_1d_ < 0) {
    /* Give debugging details if we have bad parameters and throw error. Note
       that _n_cells_1d_ can be 0 if the interaction distance is greater than
       the system diameter */
    Logger::Info("n_cells_1d: %d, min_cell_length: %2.6f, system_radius:"
                  " %2.2f",
                  _n_cells_1d_, _min_cell_length_, system_radius);
    Logger::Error("Cell list received bad initialization parameters!");
  } else if (_n_dim_ == _n_periodic_ && _n_cells_1d_ <= 3) {
    /* Warn users that we will be effectively doing all-pair interactions. If we
       have periodic BCs, 3 cells per side will result in all pairs */
    Logger::Warning("Cell list initialized with 3 cells or fewer. Calculating"
                    " all-pair interactions");
    // Fewer than three cells will have equivalent behavior
    _n_cells_1d_ = 1;
  } else if (_n_cells_1d_ <= 2) {
    /* If we do not have periodic BCs, this warning should only apply when we
       have 2 cells per side*/
    Logger::Warning("Cell list initialized with 2 cells or fewer. Calculating"
                    " all-pair interactions");
    // Fewer than 2 cells will have equivalent behavior
    _n_cells_1d_ = 2;
  }
  _cell_length_ = (double)2 * system_radius / _n_cells_1d_;
  _n_dim_ = n_dim;
  _n_periodic_ = n_periodic;
  Logger::Info("Cell list initialized with %d cells per side of length %2.2f",
                _n_cells_1d_, _cell_length_);
}

void CellList::BuildCellList() {
  Logger::Trace("Building cell list");
  AllocateCells();
  LabelCells();
  ClearCellObjects();
  /* Build with redundant neighbor pairs for fast overlap checking of new
     objects added to cell list */
  AssignCellNeighbors(true);
}

void CellList::AllocateCells() {
  Logger::Debug("Allocating cell list");
  Logger::Trace("cell_length: %2.2f", _cell_length_);
  Logger::Trace("n_cells_1d: %d", _n_cells_1d_);

  int third_dim = (_n_dim_ == 3 ? _n_cells_1d_ : 1);
  cell_ = new Cell **[_n_cells_1d_];
  for (int i = 0; i < _n_cells_1d_; ++i) {
    cell_[i] = new Cell *[_n_cells_1d_];
    for (int j = 0; j < _n_cells_1d_; ++j) {
      cell_[i][j] = new Cell[third_dim];
    }
  }
}

void CellList::Clear() {
  if (_no_init_) return;
  _n_cells_1d_ = -1;
  _no_init_ = true;
  ClearCellObjects();
  ClearCellNeighbors();
  DeallocateCells();
}

void CellList::ResetNeighbors() {
  /* Rebuild cell neighbors without redundant neighbor pairs */
  ClearCellNeighbors();
  AssignCellNeighbors(false);
}

void CellList::DeallocateCells() {
  Logger::Debug("Deallocating cell list");
  for (int i = 0; i < _n_cells_1d_; ++i) {
    for (int j = 0; j < _n_cells_1d_; ++j) {
      delete[] cell_[i][j];
    }
    delete[] cell_[i];
  }
  delete[] cell_;
}

void CellList::LabelCells() {
  int third_dim = (_n_dim_ == 3 ? _n_cells_1d_ : 1);
  for (int i = 0; i < _n_cells_1d_; ++i) {
    for (int j = 0; j < _n_cells_1d_; ++j) {
      for (int k = 0; k < third_dim; ++k) {
        cell_[i][j][k].AssignIndex(i, j, k);
      }
    }
  }
}

void CellList::MakePairs(std::vector<Interaction> &pair_list) {
  Logger::Debug("Constructing object interaction pairs");
  int third_dim = (_n_dim_ == 3 ? _n_cells_1d_ : 1);
  for (int i = 0; i < _n_cells_1d_; ++i) {
    for (int j = 0; j < _n_cells_1d_; ++j) {
      for (int k = 0; k < third_dim; ++k) {
        cell_[i][j][k].MakePairs(pair_list);
      }
    }
  }
}

xyz_coord CellList::FindCellCoords(Object &obj) {
  //const double *const spos = obj.GetScaledPosition();
  const double *const pos = obj.GetPosition();

  /*if (1 == 1) {
  double x = pos[0]/74 + 0.5;
  int xcell = (int)floor(_n_cells_1d_ * x);
  if (xcell == _n_cells_1d_)
    xcell -= 1;
  double y = pos[1]/74 + 0.5;
  int ycell = (int)floor(_n_cells_1d_ * y);
  if (ycell == _n_cells_1d_)
    ycell -= 1;
  int zcell = 0;
  if (_n_dim_ == 3) {
    double z = pos[2]/74 + 0.5;
    zcell = (int)floor(_n_cells_1d_ * z);
    if (zcell == _n_cells_1d_)
      zcell -= 1;
  }
  return std::make_tuple(xcell, ycell, zcell);
  }*/
  //else {
  double x = 0;
  if (pos[0]<=-159) { x = 0; }
  else if (pos[0]>=79) { x = 1; } 
  else {
  x = (pos[0]+160)*0.0042;
   }
  int xcell = (int)floor(_n_cells_1d_ * x);
  if (xcell == _n_cells_1d_)
    xcell -= 1;
  double y =0;
  if (pos[1]<=-79) { y = 0; }
  else if (pos[1]>=79) { y = .66; }  
  else {y = (pos[1] + 80)*0.0042;}
  int ycell = (int)floor(_n_cells_1d_ * y);
  if (ycell == _n_cells_1d_)
    ycell -= 1;
  int zcell = 0;
  if (_n_dim_ == 3) {
    double z =0;
    if (pos[2]<=-79) { z = 0; }
    else if (pos[2]>=79) { z = .66; } 
 
    else {z = (pos[2] + 80)*0.0042;}
    zcell = (int)floor(_n_cells_1d_ * z);
    if (zcell == _n_cells_1d_)
      zcell -= 1;
  }
  //Logger::Info("Cells are %i, %i, %i, with positions, %f, %f,%f", xcell, ycell, zcell, pos[0], pos[1], pos[2]);
  return std::make_tuple(xcell, ycell, zcell);
  //}
  //double x_ratio=pos[0]/spos[0];
  //double y_ratio=pos[1]/spos[1];
  //Logger::Info("Coordinates are %f, %f", x_ratio, y_ratio);
  Logger::Info("Cells are %i, %i, %i", xcell, ycell, zcell);
  //return std::make_tuple(xcell, ycell, zcell);
}

void CellList::ClearCellObjects() {
  Logger::Trace("Clearing cell list objects");
  int third_dim = (_n_dim_ == 3 ? _n_cells_1d_ : 1);
  for (int i = 0; i < _n_cells_1d_; ++i) {
    for (int j = 0; j < _n_cells_1d_; ++j) {
      for (int k = 0; k < third_dim; ++k) {
        cell_[i][j][k].ClearObjs();
      }
    }
  }
}

void CellList::ClearCellNeighbors() {
  Logger::Trace("Clearing cell list neighbors");
  int third_dim = (_n_dim_ == 3 ? _n_cells_1d_ : 1);
  for (int i = 0; i < _n_cells_1d_; ++i) {
    for (int j = 0; j < _n_cells_1d_; ++j) {
      for (int k = 0; k < third_dim; ++k) {
        cell_[i][j][k].ClearNeighbors();
      }
    }
  }
}

void CellList::RenewObjectsCells(std::vector<Object *> &objs) {
  ClearCellObjects();
  AssignObjectsCells(objs);
}

void CellList::AssignObjectsCells(std::vector<Object *> &objs) {
  Logger::Debug("Assigning objects to cells");
  for (auto obj = objs.begin(); obj != objs.end(); ++obj) {
    int x, y, z;
    std::tie(x, y, z) = FindCellCoords(**obj);
//#ifdef TRACE
    //Logger::Info("Object %d with type %s assigned to %s", (*obj)->GetOID(),((*obj)->GetName()).c_str(),
    //              cell_[x][y][z].Report().c_str());
//#endif
    cell_[x][y][z].AddObj(**obj);
  }
}

/* In order to have a cell list with unique cell neighbors, we want:
   For a given cell at position x_i, y_i, z_i, we want all 9 cells that are
   adjacent to us at (*, *, z_{i+1}), the three y cells adjacent to us at
   (*, y_{i+1}, z_i), and one cell adjacent to us along x at (x_{i+1}, y_i, z_i)
 */
void CellList::AssignCellNeighbors(bool redundancy) {
  /* If we build with redundancy, all cells will store neighbors of all
     adjacent cells, so total neighbor pairs stored will be doubled. This may
     be useful for quickly determining the potential interactions from a
     single object (ie useful for quick overlap checking of new objects) */
  if (redundancy) {
    Logger::Debug(
        "Assigning cell list neighbors with redundant neighbor pairs");
  } else {
    Logger::Debug("Assigning cell list neighbors");
  }
  int third_dim = (_n_dim_ == 3 ? _n_cells_1d_ : 1);
  // Loop through all cells in cell list
  for (int z = 0; z < third_dim; ++z) {
    for (int y = 0; y < _n_cells_1d_; ++y) {
      for (int x = 0; x < _n_cells_1d_; ++x) {
        Cell &c = cell_[x][y][z];
        int z_begin = (redundancy && _n_dim_ == 3 ? z - 1 : z);
        int z_end = (_n_dim_ == 3 ? z + 2 : z + 1);
        /* Add all adjacent cells "above" this cell along z axis */
        for (int zp = z_begin; zp < z_end; ++zp) {
          int nz = zp;
          if ((nz < 0 || nz == _n_cells_1d_) && _n_periodic_ >= 3) {
            nz = (nz < 0 ? _n_cells_1d_ - 1 : 0);
          } else if (nz < 0 || nz == _n_cells_1d_) {
            continue;
          }
          int y_begin = (redundancy || zp > z ? y - 1 : y);
          for (int yp = y_begin; yp < y + 2; ++yp) {
            int ny = yp;
            if ((ny < 0 || ny == _n_cells_1d_) && _n_periodic_ >= 2) {
              ny = (ny < 0 ? _n_cells_1d_ - 1 : 0);
            } else if (ny < 0 || ny == _n_cells_1d_) {
              continue;
            }
            int x_begin = (redundancy || yp > y || zp > z ? x - 1 : x + 1);
            for (int xp = x_begin; xp < x + 2; ++xp) {
              int nx = xp;
              if ((nx < 0 || nx == _n_cells_1d_) && _n_periodic_ >= 1) {
                nx = (nx < 0 ? _n_cells_1d_ - 1 : 0);
              } else if (nx < 0 || nx == _n_cells_1d_) {
                continue;
              }
              /* Note that cells never add themselves */
              if (nx == x && ny == y && nz == z) {
                continue;
              }
              c.AddNeighbor(cell_[nx][ny][nz]);
              Logger::Trace("%s has neighbor %s", c.Report().c_str(),
                            cell_[nx][ny][nz].Report().c_str());
            }
          }
        }
      }
    }
  }
}

/* Relies on redundant cell list neighbor pairs to identify all potential
   interactions with this object */
void CellList::PairSingleObject(Object &obj,
                                std::vector<Interaction> &pair_list) {
  int x, y, z;
  std::tie(x, y, z) = FindCellCoords(obj);
  Logger::Trace("Making pairs with single object %d in %s", obj.GetOID(),
                cell_[x][y][z].Report().c_str());
  const std::vector<Cell *> neighbors = cell_[x][y][z].GetCellNeighbors();
  cell_[x][y][z].PairSingleObject(obj, pair_list);
  for (auto cell = neighbors.begin(); cell != neighbors.end(); ++cell) {
    (*cell)->PairSingleObject(obj, pair_list);
  }
}
