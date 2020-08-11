#ifndef NOGRAPH

#include "cglass/graphics.hpp"
#include "cglass/macros.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdarg.h>

static void key_callback2d(GLFWwindow *window, int key, int scancode,
                           int action, int mods) {
  // set the window pointer to the keyboard object
  // You must set this pointer first before key callbacks are called!
  Graphics *graphics =
      static_cast<Graphics *>(glfwGetWindowUserPointer(window));

  if (action == GLFW_PRESS) {
    switch (key) {
    case GLFW_KEY_ESCAPE:
      graphics->keep_going_ = 1;
      break;
    case GLFW_KEY_UP:
      graphics->yTrans_ += 0.5;
      break;
    case GLFW_KEY_DOWN:
      graphics->yTrans_ -= 0.5;
      break;
    case GLFW_KEY_RIGHT:
      graphics->xTrans_ += 0.5;
      break;
    case GLFW_KEY_LEFT:
      graphics->xTrans_ -= 0.5;
      break;
    case GLFW_KEY_Z:
      if (mods & GLFW_MOD_SHIFT)
        graphics->zAngle_ += 5.0;
      else
        graphics->zAngle_ -= 5.0;
      break;
    case GLFW_KEY_SPACE:
      graphics->xyzScale_ += 1.0 / 10.0;
      break;
    case GLFW_KEY_BACKSPACE:
      graphics->xyzScale_ -= 1.0 / 10.0;
      break;
    case GLFW_KEY_C:
      if (mods & GLFW_MOD_SHIFT)
        graphics->color_switch_ = (graphics->color_switch_ + 1) % 3;
      else {
        graphics->color_switch_ = (graphics->color_switch_ - 1);
        if (graphics->color_switch_ < 0)
          graphics->color_switch_ = 2;
      }
      break;
    case GLFW_KEY_A:
      if (mods & GLFW_MOD_SHIFT) {
        graphics->alpha_ += 0.1;
        if (graphics->alpha_ > 1.0)
          graphics->alpha_ = 1.0;
      } else {
        graphics->alpha_ -= 0.1;
        if (graphics->alpha_ < 0.0)
          graphics->alpha_ = 0.0;
      }
      break;
    case GLFW_KEY_T:
      // graphics->draw_trajectory_flag_ =
      //   (graphics->draw_trajectory_flag_ + 1) % 2;
      break;
    default:
      break;
    }
  }
}

static void key_callback3d(GLFWwindow *window, int key, int scancode,
                           int action, int mods) {
  // set the window pointer to the keyboard object
  // You must set this pointer first before key callbacks are called!
  Graphics *graphics =
      static_cast<Graphics *>(glfwGetWindowUserPointer(window));

  if (action == GLFW_PRESS) {
    switch (key) {
    case GLFW_KEY_ESCAPE:
      graphics->keep_going_ = 1;
      break;
    case GLFW_KEY_UP:
      graphics->xAngle_ += 5.0;
      break;
    case GLFW_KEY_DOWN:
      graphics->xAngle_ -= 5.0;
      break;
    case GLFW_KEY_RIGHT:
      graphics->yAngle_ += 5.0;
      break;
    case GLFW_KEY_LEFT:
      graphics->yAngle_ -= 5.0;
      break;
    case GLFW_KEY_Z:
      if (mods & GLFW_MOD_SHIFT)
        graphics->zAngle_ -= 5.0;
      else
        graphics->zAngle_ += 5.0;
      break;
    case GLFW_KEY_C:
      if (mods & GLFW_MOD_SHIFT)
        graphics->cNear_ -= 0.5;
      else
        graphics->cNear_ += 0.5;
      break;
    case GLFW_KEY_V:
      if (mods & GLFW_MOD_SHIFT)
        graphics->cFar_ += 0.5;
      else
        graphics->cFar_ -= 0.5;
      break;
    case GLFW_KEY_J:
      if (mods & GLFW_MOD_SHIFT)
        graphics->cxAngle_ -= 5.0;
      else
        graphics->cxAngle_ += 5.0;
      break;
    case GLFW_KEY_K:
      if (mods & GLFW_MOD_SHIFT)
        graphics->cyAngle_ -= 5.0;
      else
        graphics->cyAngle_ += 5.0;
      break;
    case GLFW_KEY_SPACE:
      graphics->xyzScale_ += 1.0 / 10.0;
      break;
    case GLFW_KEY_BACKSPACE:
      graphics->xyzScale_ -= 1.0 / 10.0;
      break;
    case GLFW_KEY_M:
      graphics->color_switch_ = (graphics->color_switch_ + 1) % 3;
      break;
    case GLFW_KEY_R:
      graphics->xyzScale_ *= 0.5;
      break;
    case GLFW_KEY_T:
      graphics->xyzScale_ *= 2.0;
      break;
    default:
      break;
    }
  }
}

void Graphics::KeyInteraction() {
  // Tell glfw where the graphics class is
  glfwSetWindowUserPointer(window_, static_cast<void *>(this));
  // Set the proper key callback, since 2d and 3d should have different key
  // controls
  // No key interactions if we're doing auto_graph
  //if (!auto_graph_) {
  if (n_dim_ == 2)
    glfwSetKeyCallback(window_, key_callback2d);
  else if (n_dim_ == 3)
    glfwSetKeyCallback(window_, key_callback3d);
  //}
  // Check window for key presses/mouse interaction
  glfwPollEvents();
  // Exit program if you close graphics window
  if (glfwWindowShouldClose(window_)) {
    glfwDestroyWindow(window_);
    glfwTerminate();
    Logger::Error("Graphics window closed: terminating program");
  }
}

void Graphics::ScalePositions() {
  for (auto it = graph_array_->begin(); it != graph_array_->end(); ++it)
    for (int i = 0; i < space_->n_periodic; ++i)
      (*it)->r[i] = unit_cell_[n_dim_ * i + i] * (*it)->r[i];
}

void Graphics::Init(std::vector<graph_struct *> *graph_array,
                    space_struct *s_struct, double background,
                    int draw_boundary, int auto_graph, double alpha) {
  space_ = s_struct;
  graph_array_ = graph_array;
  n_dim_ = space_->n_dim;
  unit_cell_ = space_->unit_cell;
  auto_graph_ = auto_graph;
  alpha_ = alpha;
  // ScalePositions();
  if (draw_boundary == 0 || space_->type == +boundary_type::wall) {
    boundary_ = boundary_type::none;
  } else {
    boundary_ = space_->type;
  }
  for (int i = 0; i < 3; ++i) {
    background_color_[i] = background;
  }
  z_correct_ = 0;
  if (boundary_ == +boundary_type::budding) {
    double m_rad = space_->radius;
    double d_rad = space_->bud_radius;
    double m_d_dist = space_->bud_height;
    z_correct_ = 0.5 * (m_d_dist + m_rad + d_rad) - m_rad;
  }
  keep_going_ = 0;
  InitColormap();
  InitWindow();
  InitDiscoRectangle();
  InitSpheroCylinder();
  // InitText();
}

void Graphics::Clear() {
  keep_going_ = 0;
  delete[] discorectangle_.vertex_buffer_data_;
  delete[] discorectangle_.element_buffer_data_;
  delete[] spherocylinder_.vertex_buffer_data_;
  delete[] spherocylinder_.element_buffer_data_;
  glfwDestroyWindow(window_);
  glfwTerminate();
}

void Graphics::InitWindow() {
  if (n_dim_ == 2)
    Init2dWindow();
  else if (n_dim_ == 3)
    Init3dWindow();
}

void Graphics::Init2dWindow() {
  // Defaults for window view and coloring
  xyzScale_ = 0.95;
  xTrans_ = yTrans_ = 0.0; // x,y translation
  color_switch_ = 1;       // Default 'nice' spherocylinder coloring
  //alpha_ = 1.0;            // transparency

  /* Compute size of the window respecting the aspect ratio. */
  if (unit_cell_[0] > unit_cell_[3]) {
    windx_ = 800;
    windy_ = (int)(windx_ * unit_cell_[3] / unit_cell_[0]);
  } else {
    windy_ = 800;
    windx_ = (int)(windy_ * unit_cell_[0] / unit_cell_[3]);
  }

  { // Make dummy window so we can get GL extensions
    glfwWindowHint(GLFW_VISIBLE, GL_FALSE);
    if (!glfwInit())
      Logger::Error("Unable to initialize GLFW toolkit\n");

    window_ = glfwCreateWindow(windx_, windy_, "2D Sphero", NULL, NULL);

    if (!window_) {
      glfwTerminate();
      Logger::Error("Failed to create display window\n");
    }

    // Wrangle in GL function pointers for GL version > 2.0
    glewInit();
    glfwDestroyWindow(window_);
  }

  /* Rebuild window now that we have proper information about available GL
   * extensions */
  glfwWindowHint(GLFW_SAMPLES, 16);
  window_ = glfwCreateWindow(windx_, windy_, "2D Sphero", NULL, NULL);

  if (!window_) {
    glfwTerminate();
    Logger::Error("Failed to create display window\n");
  }

  glfwMakeContextCurrent(window_);
  glfwSwapInterval(0); // swap buffers immediately
  glewInit();          // re-initialize glew for current window context

  // Blending and antialiasing. Makes things pretty.
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  glEnable(GL_POLYGON_SMOOTH);
  glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
  glEnable(GL_MULTISAMPLE);

  /* frame buffer clears should be to white by default in 2d */
  glClearColor(background_color_[0], background_color_[1], background_color_[2],
               1.0);

  /* set up projection transform */
  double xmin = -unit_cell_[0] / 2.;
  double xmax = unit_cell_[0] / 2.;
  double ymin = -unit_cell_[3] / 2.;
  double ymax = unit_cell_[3] / 2.;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(
      (GLdouble)(xmin - 0.05 * unit_cell_[0]),
      (GLdouble)(xmax + 0.05 * unit_cell_[0]),
      (GLdouble)(ymin - 0.05 * unit_cell_[3]),
      (GLdouble)(ymax +
                 0.05 *
                     unit_cell_[3])); // deprecated, but annoying to implement
  /* JMM attempting to fix Mac resize issues */
  glfwGetFramebufferSize(window_, &windx_, &windy_);
  /* establish initial viewport */
  glViewport(0, 0, windx_, windy_);
}

void Graphics::Init3dWindow() {
  xAngle_ = 90;
  yAngle_ = 0;
  zAngle_ = 90;
  xyzScale_ = 0.95; // zoom out just a hair.
  zTrans_ = -24;    // camera position in z
  color_switch_ = 1;

  // FIXME: Assumes a square box. Better way to calculate this, but we
  // typically use square boxes for now
  double a_perp_max = 0.0;
  for (int i = 0; i < 3; ++i)
    a_perp_max = MAX(unit_cell_[3 * i + i], a_perp_max);

  cNear_ = -2.0 * a_perp_max; // clipping plane near
  cFar_ = 2.0 * a_perp_max;   // clipping plane far
  cxAngle_ = cyAngle_ = 0;    // plane angles

  windx_ = windy_ = 400; // window size in pixels

  { // Make dummy window so we can get GL extensions
    glfwWindowHint(GLFW_VISIBLE, GL_FALSE);
    if (!glfwInit())
      Logger::Error("Unable to initialize GLFW toolkit\n");

    window_ = glfwCreateWindow(windx_, windy_, "2D Sphero", NULL, NULL);

    if (!window_) {
      glfwTerminate();
      Logger::Error("Failed to create display window\n");
    }
    /* JMM - attempting to fix window size issue with mac Retina display */
    glfwGetFramebufferSize(window_, &windx_, &windy_);
    glViewport(0, 0, windx_, windy_);
    /* Wrangle in GL function pointers for GL version > 2.0 */
    glewInit();
    glfwDestroyWindow(window_);
  }

  // Rebuild window now that we have proper information about the GL extensions
  glfwWindowHint(GLFW_SAMPLES, 16); // 16x Multisampling
  window_ = glfwCreateWindow(windx_, windy_, "2D Sphero", NULL, NULL);

  if (!window_) {
    glfwTerminate();
    Logger::Error("Failed to create display window\n");
  }

  glfwMakeContextCurrent(window_);
  glfwSwapInterval(0); // Swap buffers immediately
  glewInit();          // Get GL extensions for current context

  /* FIXME: this probably does nothing nowadays */
  /* Material properties. */
  GLfloat no_mat[] = {0.0, 0.0, 0.0, 1.0};
  GLfloat mat_diffuse[] = {1.0, 1.0, 1.0, 1.0};
  GLfloat low_shininess[] = {5.0};
  glMaterialfv(GL_FRONT, GL_AMBIENT, no_mat);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
  glMaterialfv(GL_FRONT, GL_SHININESS, low_shininess);
  glMaterialfv(GL_FRONT, GL_EMISSION, no_mat);

  /* Lights. */
  GLfloat ambient[] = {0.0, 0.0, 0.0, 1.0};
  GLfloat diffuse[] = {1.0, 1.0, 1.0, 1.0};
  GLfloat position[] = {0.0, 2.0, 3.0, 0.0};
  GLfloat lmodel_ambient[] = {0.4, 0.4, 0.4, 1.0};
  GLfloat local_view[] = {0.0};
  glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
  glLightfv(GL_LIGHT0, GL_POSITION, position);
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
  glLightModelfv(GL_LIGHT_MODEL_LOCAL_VIEWER, local_view);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  /* Used in order to change one material property (here the color). */
  glColorMaterial(GL_FRONT, GL_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);

  /* Blending and antialiasing. */
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
  glEnable(GL_MULTISAMPLE);

  zTrans_ = -1.5 * unit_cell_[8];
  glEnable(GL_DEPTH_TEST); /* enable depth buffering */
  glDepthFunc(GL_LESS);    /* pedantic, GL_LESS is the default */
  glClearDepth(1.0);       /* pedantic, 1.0 is the default */

  /* frame buffer clears should be to white */
  glClearColor(background_color_[0], background_color_[1], background_color_[2],
               1); /* FIXME: Guess that we're drawing
         spherical confined system */
  /* set up projection transform */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(-60.0, 1.0, 1.5 * unit_cell_[8] / 2.0,
                 10.0 * 1.5 * unit_cell_[8] / 2.0); // deprecated

  /* establish initial viewport */
  glfwGetFramebufferSize(window_, &windx_, &windy_);
  glViewport(0, 0, windx_, windy_);
}

void Graphics::InitColormap() { // Pretty straight forward. Builds
                                // colormap. Automatic for all instances of
                                // Graphics class after Init() call
  /* Loop over colormap indices. */
  for (int i_color = 0; i_color < n_rgb_; ++i_color) {
    GLfloat red, green, blue;
    if (i_color < 64) {
      red = 1.0;
      green = i_color / 64.0;
      blue = 0.0;
    } else if (i_color >= 64 && i_color < 128) {
      red = (128.0 - i_color) / 64.0;
      green = 1.0;
      blue = 0.0;
    } else if (i_color >= 128 && i_color < 192) {
      red = 0.0;
      green = 1.0;
      blue = (i_color - 128.0) / 64.0;
    } else if (i_color >= 192 && i_color < 256) {
      red = 0.0;
      green = (256.0 - i_color) / 64.0;
      blue = 1.0;
    } else if (i_color >= 256 && i_color < 320) {
      red = (i_color - 256.0) / 64.0;
      green = 0.0;
      blue = 1.0;
    } else {
      red = 1.0;
      green = 0.0;
      blue = (384.0 - i_color) / 64.0;
    }

    colormap_[i_color * 4 + 0] = red;
    colormap_[i_color * 4 + 1] = green;
    colormap_[i_color * 4 + 2] = blue;
    colormap_[i_color * 4 + 3] = 1.0;
  }
}

void Graphics::DrawDiscorectangles() {
  glMatrixMode(GL_MODELVIEW); // Make sure we're using the model transform

  // Tells gl to use the appropriate shader for the discorectangle
  glUseProgram(discorectangle_.program_); // This could be moved to a
                                          // GraphicsPrimitive method.

  // spherocylinder takes two parameters, half its length, and its diameter
  // Get location of half_l and diameter parameters, this shouldn't need to be
  // done everytime but for whatever reason it wasn't grabbing the correct
  // location at initialization
  discorectangle_.uniforms_.half_l =
      glGetUniformLocation(discorectangle_.program_, "half_l");
  discorectangle_.uniforms_.diameter =
      glGetUniformLocation(discorectangle_.program_, "diameter");

  // Prep to send vertex data to shader
  glEnableVertexAttribArray(discorectangle_.attributes_.position);
  glBindBuffer(GL_ARRAY_BUFFER, discorectangle_.vertex_buffer_);
  glVertexAttribPointer(discorectangle_.attributes_.position, // attribute
                        2,        // number of elements per vertex, here (x,y,z)
                        GL_FLOAT, // the type of each element
                        GL_FALSE, // take our values as-is
                        0,        // no extra data between each position
                        (void *)0 // offset of first element
  );

  // Use the element buffer (numerical pointer to each set of vertices that make
  // a triangle)
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, discorectangle_.element_buffer_);

  // Set default bond color
  GLfloat color[4] = {0.0, 0.0, 1.0, 1.0};
  for (auto it = graph_array_->begin(); it != graph_array_->end(); ++it) {
    double theta = atan2((*it)->u[1], (*it)->u[0]); // rotation angle
    double theta_color = 0;
    // Check draw type
    if ((*it)->draw == +draw_type::fixed) {
      // draw is set to flat
      theta_color = (*it)->color;
    } else if ((*it)->draw == +draw_type::orientation) {
      // orientation draw
      theta_color = theta;
    }
    if ((*it)->draw == +draw_type::bw) {
      for (int i = 0; i < 3; ++i)
        color[i] = (background_color_[0] < 1 ? 1 : 0.1);
      glColor4fv(color);
    } else {
      theta_color = theta_color + 1.5 * M_PI;
      while (theta_color > 2.0 * M_PI) {
        theta_color -= 2.0 * M_PI;
      }

      // Color based on orientation
      int color_index = (int)(((theta_color) / (2.0 * M_PI)) * n_rgb_);
      if (color_index < 0)
        color_index = 0;
      else if (color_index >= n_rgb_)
        color_index = n_rgb_ - 1;

      glColor4fv(&colormap_[color_index * 4]);
    }
    glPushMatrix(); // Duplicate current modelview
    {
      glTranslatef((*it)->r[0], (*it)->r[1] - z_correct_,
                   0.0); // Translate rod
      glRotatef((GLfloat)((theta / M_PI) * 180.0 - 90.0), 0.0, 0.0,
                1.0); // rotate rod

      // Tell shader about our spherocylinder parameters
      double half_length = 0.5 * (*it)->length;
      glUniform1f(discorectangle_.uniforms_.half_l, half_length);
      glUniform1f(discorectangle_.uniforms_.diameter, (*it)->diameter);
      glDrawElements(GL_TRIANGLES, discorectangle_.n_triangles_ * 3,
                     GL_UNSIGNED_SHORT,
                     (void *)0); // draw.
    }
    glPopMatrix(); // restore default modelview
  }
  glDisableVertexAttribArray(
      discorectangle_.attributes_
          .position); // Tell GL to forget about our vertex array
}

void Graphics::Draw2d() {
  KeyInteraction();
  UpdateWindow(); // Recalculate window parameters (in case of resize)
  DrawDiscorectangles();
  DrawBoundary();
  // DrawText();
  glfwSwapBuffers(window_); // Print frame to screen
}

// void Graphics::DrawText() {

// glMatrixMode(GL_MODELVIEW); // Make sure we're using the model transform
// glUseProgram(text_.program_);

//// In case things are not working right
// text_.attribute_coord_ = glGetAttribLocation(text_.program_, "coord");
// text_.uniform_tex_ = glGetUniformLocation(text_.program_, "tex");
// text_.uniform_color_ = glGetUniformLocation(text_.program_, "color");

// glVertexAttribPointer(text_.attribute_coord_, // attribute
// 2,                 // number of elements per vertex, here (x,y,z)
// GL_FLOAT,          // the type of each element
// GL_FALSE,          // take our values as-is
// 0,                 // no extra data between each position
//(void*)0                  // offset of first element
//);

//// Use the element buffer (numerical pointer to each set of vertices that make
/// a triangle)
////glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, text_.element_buffer_);

// GLfloat black[4] = { 0, 0, 0, 1 };
// GLfloat red[4] = { 1, 0, 0, 1 };
// GLfloat transparent_green[4] = { 0, 1, 0, 0.5 };

//[> Set font size to 48 pixels, color to black <]
// FT_Set_Pixel_Sizes(text_.face_, 0, 48);
// glUniform4fv(text_.uniform_color_, 1, black);

//[> Effects of alignment <]
// text_.RenderText("The Quick Brown Fox Jumps Over The Lazy Dog", 0,0, 0, 0);
// text_.RenderText("The Quick Brown Fox Jumps Over The Lazy Dog", -1 + 8 *
// windx_, 1 - 50 * windy_, windx_, windy_); text_.RenderText("The Misaligned
// Fox Jumps Over The Lazy Dog", -1 + 8.5 * windx_, 1 - 100.5 * windy_, windx_,
// windy_);

//[> Scaling the texture versus changing the font size <]
// text_.RenderText("The Small Texture Scaled Fox Jumps Over The Lazy Dog", -1 +
// 8 * windx_, 1 - 175 * windy_, windx_ * 0.5, windy_ * 0.5);
// FT_Set_Pixel_Sizes(text_.face_, 0, 24);
// text_.RenderText("The Small Font Sized Fox Jumps Over The Lazy Dog", -1 + 8 *
// windx_, 1 - 200 * windy_, windx_, windy_); FT_Set_Pixel_Sizes(text_.face_, 0,
// 48); text_.RenderText("The Tiny Texture Scaled Fox Jumps Over The Lazy Dog",
// -1 + 8 * windx_, 1 - 235 * windy_, windx_ * 0.25, windy_ * 0.25);
// FT_Set_Pixel_Sizes(text_.face_, 0, 12);
// text_.RenderText("The Tiny Font Sized Fox Jumps Over The Lazy Dog", -1 + 8 *
// windx_, 1 - 250 * windy_, windx_, windy_); FT_Set_Pixel_Sizes(text_.face_, 0,
// 48);

//[> Colors and transparency <]
// text_.RenderText("The Solid Black Fox Jumps Over The Lazy Dog", -1 + 8 *
// windx_, 1 - 430 * windy_, windx_, windy_);

// glUniform4fv(text_.uniform_color_, 1, red);
// text_.RenderText("The Solid Red Fox Jumps Over The Lazy Dog", -1 + 8 *
// windx_, 1 - 330 * windy_, windx_, windy_); text_.RenderText("The Solid Red
// Fox Jumps Over The Lazy Dog", -1 + 28 * windx_, 1 - 450 * windy_, windx_,
// windy_);

// glUniform4fv(text_.uniform_color_, 1, transparent_green);
// text_.RenderText("The Transparent Green Fox Jumps Over The Lazy Dog", -1 + 8
// * windx_, 1 - 380 * windy_, windx_, windy_); text_.RenderText("The
// Transparent Green Fox Jumps Over The Lazy Dog", -1 + 18 * windx_, 1 - 440 *
// windy_, windx_, windy_);
//}

void Graphics::DrawSpheros() {
  glMatrixMode(GL_MODELVIEW); // Use modelview matrix

  /* Don't draw back facing triangles (as they would be occluded anyway */
  glEnable(GL_CULL_FACE);
  /* Use vertex/fragment custom shader for spherocylinder */
  glUseProgram(spherocylinder_.program_);
  /* Get location of half_l parameter, this shouldn't need to be done everytime
    but for whatever reason it wasn't grabbing the correct location at
    initialization */
  spherocylinder_.uniforms_.half_l =
      glGetUniformLocation(spherocylinder_.program_, "half_l");
  spherocylinder_.uniforms_.diameter =
      glGetUniformLocation(spherocylinder_.program_, "diameter");
  /* Prep to send vertex data to shader */
  glEnableVertexAttribArray(spherocylinder_.attributes_.position);
  glBindBuffer(GL_ARRAY_BUFFER, spherocylinder_.vertex_buffer_);
  glVertexAttribPointer(spherocylinder_.attributes_.position, // attribute
                        3,        // number of elements per vertex, here (x,y,z)
                        GL_FLOAT, // the type of each element
                        GL_FALSE, // take our values as-is
                        0,        // no extra data between each position
                        (void *)0 // offset of first element
  );

  /* Use the element buffer (numerical pointer to each set of vertices that make
   * a triangle) */
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, spherocylinder_.element_buffer_);

  GLfloat color[4] = {0.0, 0.0, 1.0, 1.0}; // default bond color
  // for (int i_bond = 0; i_bond < n_spheros; ++i_bond) {
  for (auto it = graph_array_->begin(); it != graph_array_->end(); ++it) {
    /* Determine phi rotation angle, amount to rotate about y. */
    double phi = acos((*it)->u[2]);
    double theta = 0.0;
    double theta_color = 0.0;
    /* Determine theta rotation, amount to rotate about z. */
    double length_xy = sqrt(SQR((*it)->u[0]) + SQR((*it)->u[1]));
    if (length_xy > 0.0) {
      theta = acos((*it)->u[0] / length_xy);
      if ((*it)->u[1] < 0.0)
        theta = (2.0 * M_PI) - theta;
    }
    color[3] = alpha_;
    std::string d_type((*it)->draw._to_string());
    // Color is now based on draw for the object
    if (d_type.compare("fixed") == 0) {
      // flat color
      theta_color = (*it)->color;
    } else if (d_type.compare("orientation") == 0) {
      theta_color = theta;
    }

    if (d_type.compare("bw") == 0) {
      for (int i = 0; i < 3; ++i)
        color[i] = (background_color_[0] < 1 ? 1 : 0.1);
      glColor4fv(color);
    } else {
      /* Convert from HSL to RGB coloring scheme, unique orientation coloring
       * scheme */
      double L = 0.3 * (*it)->u[2] + 0.5;
      double C = (1 - ABS(2 * L - 1));
      double H_prime = 3.0 * theta_color / M_PI;
      double X = C * (1.0 - ABS(fmod(H_prime, 2.0) - 1));

      if (H_prime < 0.0) {
        color[0] = 0.0;
        color[1] = 0.0;
        color[2] = 0.0;
      } else if (H_prime < 1.0) {
        color[0] = C;
        color[1] = X;
        color[2] = 0;
      } else if (H_prime < 2.0) {
        color[0] = X;
        color[1] = C;
        color[2] = 0;
      } else if (H_prime < 3.0) {
        color[0] = 0;
        color[1] = C;
        color[2] = X;
      } else if (H_prime < 4.0) {
        color[0] = 0;
        color[1] = X;
        color[2] = C;
      } else if (H_prime < 5.0) {
        color[0] = X;
        color[1] = 0;
        color[2] = C;
      } else if (H_prime < 6.0) {
        color[0] = C;
        color[1] = 0;
        color[2] = X;
      }

      double m = L - 0.5 * C;
      color[0] = color[0] + m;
      color[1] = color[1] + m;
      color[2] = color[2] + m;

      glColor4fv(color);
    }
    /* Get position of spherocylinder center. */
    GLfloat v0 = (*it)->r[0];
    GLfloat v1 = (*it)->r[1];
    GLfloat v2 = (*it)->r[2] - z_correct_;

    /* Let the shader know the length of our bond */
    GLfloat half_length = 0.5 * (*it)->length;
    glUniform1f(spherocylinder_.uniforms_.half_l, half_length);
    glUniform1f(spherocylinder_.uniforms_.diameter, (*it)->diameter);
    /* Make copy of modelview matrix to work on */
    glPushMatrix();

    /* Move and rotate originally z oriented, origin centered bond to new angle
     * and position */
    glTranslatef(v0, v1, v2);
    if (theta != 0.0)
      glRotatef((GLfloat)((theta / M_PI) * 180.0), 0.0, 0.0, 1.0);
    if (phi != 0.0)
      glRotatef((GLfloat)((phi / M_PI) * 180.0), 0.0, 1.0, 0.0);

    /* Send data to shader to be drawn */
    glDrawElements(GL_TRIANGLES, spherocylinder_.n_triangles_ * 3,
                   GL_UNSIGNED_SHORT, (void *)0);
    /* Reset modelview matrix */
    glPopMatrix();
  }
  glDisableVertexAttribArray(spherocylinder_.attributes_.position);
}

void Graphics::Draw3d() {
  KeyInteraction();
  UpdateWindow();
  DrawSpheros();
  DrawBoundary();
  // DrawText();
  glfwSwapBuffers(window_);
}

void Graphics::DrawBoundary() {
  if (boundary_ == +boundary_type::none)
    return;
  glUseProgram(0);

  GLfloat color[4] = {0.5, 0.5, 0.5, 1.0};
  glColor4fv(color);
  if (boundary_ == +boundary_type::sphere) {
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);

    /* Fixme: I'm not drawing anything */
    // Turn on wireframe mode
    if (n_dim_ == 2)
      // DrawWireSphere(0.5*unit_cell_[0], 3, 32);
      DrawWireSphere(0.5 * unit_cell_[0], 32, 32);
    if (n_dim_ == 3)
      DrawWireSphere(0.5 * unit_cell_[0], 16, 16);

  } else if (boundary_ == +boundary_type::box)
    DrawBox();
  else if (boundary_ == +boundary_type::budding)
    DrawBudding();
}

void Graphics::DrawWireSphere(double r, int lats, int longs) {
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glLineWidth(1.0);

  for (int i = 0; i <= lats; i++) {
    double lat0 = M_PI * (-0.5 + (double)(i - 1) / lats);
    double z0 = r * sin(lat0);
    double zr0 = cos(lat0);
    double lat1 = M_PI * (-0.5 + (double)i / lats);
    double z1 = r * sin(lat1);
    double zr1 = cos(lat1);
    glBegin(GL_QUAD_STRIP);
    for (int j = 0; j <= longs; j++) {
      double lng = 2 * M_PI * (double)(j - 1) / longs;
      double x = r * cos(lng);
      double y = r * sin(lng);
      glNormal3f(x * zr0, y * zr0, z0);
      glVertex3f(x * zr0, y * zr0, z0);
      glNormal3f(x * zr1, y * zr1, z1);
      glVertex3f(x * zr1, y * zr1, z1);
    }
    glEnd();
  }
  // Turn off wireframe mode
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glEnable(GL_LIGHTING);
  glEnable(GL_CULL_FACE);
}

void Graphics::DrawLoop() {
  if (auto_graph_) {
    Draw();
    return;
  }
  ScalePositions();
  do {
    if (n_dim_ == 2)
      Draw2d();
    else if (n_dim_ == 3)
      Draw3d();
  } while (!keep_going_);
  return;
}

void Graphics::Draw() {
  ScalePositions();
  if (n_dim_ == 2)
    Draw2d();
  else if (n_dim_ == 3)
    Draw3d();
}

void Graphics::DrawSquare() {
  glUseProgram(0);

  /* Draw the bounding box. */
  GLfloat box_color[4] = {1.0, 0.0, 0.5, 1.0};
  GLfloat v0, v1;
  GLfloat w = 2;
  glLineWidth(w);
  glColor4fv(box_color);
  glBegin(GL_LINE_LOOP);
  v0 = -unit_cell_[0] / 2 - unit_cell_[1] / 2;
  v1 = -unit_cell_[3] / 2 - unit_cell_[2] / 2;
  glVertex2f(v0, v1);
  v0 = -unit_cell_[0] / 2 + unit_cell_[1] / 2;
  v1 = unit_cell_[3] / 2 - unit_cell_[2] / 2;
  glVertex2f(v0, v1);
  v0 = unit_cell_[0] / 2 + unit_cell_[1] / 2;
  v1 = unit_cell_[3] / 2 + unit_cell_[2] / 2;
  glVertex2f(v0, v1);
  v0 = unit_cell_[0] / 2 - unit_cell_[1] / 2;
  v1 = -unit_cell_[3] / 2 + unit_cell_[2] / 2;
  glVertex2f(v0, v1);
  v0 = -unit_cell_[0] / 2 - unit_cell_[1] / 2;
  v1 = -unit_cell_[3] / 2 - unit_cell_[2] / 2;
  glVertex2f(v0, v1);
  glEnd();
}

void Graphics::DrawCube() {
  int w;
  GLfloat v0, v1, v2;

  /* Turn off the light. */
  glDisable(GL_LIGHTING);

  /* Line thickness. */
  w = 3;
  glLineWidth((GLfloat)w);

  /* Draw the bounding box. */
  /* front face */
  glBegin(GL_LINE_LOOP);
  v0 = 0.5 * (-unit_cell_[0] - unit_cell_[1] - unit_cell_[2]);
  v1 = 0.5 * (-unit_cell_[3] - unit_cell_[4] - unit_cell_[5]);
  v2 = 0.5 * (-unit_cell_[6] - unit_cell_[7] - unit_cell_[8]);
  glVertex3f(v0, v1, v2);
  v0 = 0.5 * (unit_cell_[0] - unit_cell_[1] - unit_cell_[2]);
  v1 = 0.5 * (unit_cell_[3] - unit_cell_[4] - unit_cell_[5]);
  v2 = 0.5 * (unit_cell_[6] - unit_cell_[7] - unit_cell_[8]);
  glVertex3f(v0, v1, v2);
  v0 = 0.5 * (unit_cell_[0] + unit_cell_[1] - unit_cell_[2]);
  v1 = 0.5 * (unit_cell_[3] + unit_cell_[4] - unit_cell_[5]);
  v2 = 0.5 * (unit_cell_[6] + unit_cell_[7] - unit_cell_[8]);
  glVertex3f(v0, v1, v2);
  v0 = 0.5 * (-unit_cell_[0] + unit_cell_[1] - unit_cell_[2]);
  v1 = 0.5 * (-unit_cell_[3] + unit_cell_[4] - unit_cell_[5]);
  v2 = 0.5 * (-unit_cell_[6] + unit_cell_[7] - unit_cell_[8]);
  glVertex3f(v0, v1, v2);
  v0 = 0.5 * (-unit_cell_[0] - unit_cell_[1] - unit_cell_[2]);
  v1 = 0.5 * (-unit_cell_[3] - unit_cell_[4] - unit_cell_[5]);
  v2 = 0.5 * (-unit_cell_[6] - unit_cell_[7] - unit_cell_[8]);
  glVertex3f(v0, v1, v2);
  v0 = 0.5 * (-unit_cell_[0] - unit_cell_[1] + unit_cell_[2]);
  v1 = 0.5 * (-unit_cell_[3] - unit_cell_[4] + unit_cell_[5]);
  v2 = 0.5 * (-unit_cell_[6] - unit_cell_[7] + unit_cell_[8]);
  glVertex3f(v0, v1, v2);
  v0 = 0.5 * (unit_cell_[0] - unit_cell_[1] + unit_cell_[2]);
  v1 = 0.5 * (unit_cell_[3] - unit_cell_[4] + unit_cell_[5]);
  v2 = 0.5 * (unit_cell_[6] - unit_cell_[7] + unit_cell_[8]);
  glVertex3f(v0, v1, v2);
  v0 = 0.5 * (unit_cell_[0] - unit_cell_[1] - unit_cell_[2]);
  v1 = 0.5 * (unit_cell_[3] - unit_cell_[4] - unit_cell_[5]);
  v2 = 0.5 * (unit_cell_[6] - unit_cell_[7] - unit_cell_[8]);
  glVertex3f(v0, v1, v2);
  glEnd();
  glBegin(GL_LINE_LOOP);
  v0 = 0.5 * (-unit_cell_[0] + unit_cell_[1] - unit_cell_[2]);
  v1 = 0.5 * (-unit_cell_[3] + unit_cell_[4] - unit_cell_[5]);
  v2 = 0.5 * (-unit_cell_[6] + unit_cell_[7] - unit_cell_[8]);
  glVertex3f(v0, v1, v2);
  v0 = 0.5 * (-unit_cell_[0] + unit_cell_[1] + unit_cell_[2]);
  v1 = 0.5 * (-unit_cell_[3] + unit_cell_[4] + unit_cell_[5]);
  v2 = 0.5 * (-unit_cell_[6] + unit_cell_[7] + unit_cell_[8]);
  glVertex3f(v0, v1, v2);
  v0 = 0.5 * (unit_cell_[0] + unit_cell_[1] + unit_cell_[2]);
  v1 = 0.5 * (unit_cell_[3] + unit_cell_[4] + unit_cell_[5]);
  v2 = 0.5 * (unit_cell_[6] + unit_cell_[7] + unit_cell_[8]);
  glVertex3f(v0, v1, v2);
  v0 = 0.5 * (unit_cell_[0] + unit_cell_[1] - unit_cell_[2]);
  v1 = 0.5 * (unit_cell_[3] + unit_cell_[4] - unit_cell_[5]);
  v2 = 0.5 * (unit_cell_[6] + unit_cell_[7] - unit_cell_[8]);
  glVertex3f(v0, v1, v2);
  glEnd();
  glBegin(GL_LINE_LOOP);
  v0 = 0.5 * (-unit_cell_[0] - unit_cell_[1] + unit_cell_[2]);
  v1 = 0.5 * (-unit_cell_[3] - unit_cell_[4] + unit_cell_[5]);
  v2 = 0.5 * (-unit_cell_[6] - unit_cell_[7] + unit_cell_[8]);
  glVertex3f(v0, v1, v2);
  v0 = 0.5 * (-unit_cell_[0] + unit_cell_[1] + unit_cell_[2]);
  v1 = 0.5 * (-unit_cell_[3] + unit_cell_[4] + unit_cell_[5]);
  v2 = 0.5 * (-unit_cell_[6] + unit_cell_[7] + unit_cell_[8]);
  glVertex3f(v0, v1, v2);
  glEnd();
  glBegin(GL_LINE_LOOP);
  v0 = 0.5 * (unit_cell_[0] - unit_cell_[1] + unit_cell_[2]);
  v1 = 0.5 * (unit_cell_[3] - unit_cell_[4] + unit_cell_[5]);
  v2 = 0.5 * (unit_cell_[6] - unit_cell_[7] + unit_cell_[8]);
  glVertex3f(v0, v1, v2);
  v0 = 0.5 * (unit_cell_[0] + unit_cell_[1] + unit_cell_[2]);
  v1 = 0.5 * (unit_cell_[3] + unit_cell_[4] + unit_cell_[5]);
  v2 = 0.5 * (unit_cell_[6] + unit_cell_[7] + unit_cell_[8]);
  glVertex3f(v0, v1, v2);
  glEnd();

  /* Turn on light. */
  glEnable(GL_LIGHTING);
}

void Graphics::DrawBox() {
  if (n_dim_ == 2)
    DrawSquare();
  else if (n_dim_ == 3)
    DrawCube();
}

void Graphics::DrawBudding() {
  if (n_dim_ == 2)
    Draw2dBudding();
  else if (n_dim_ == 3)
    Draw3dBudding();
}

void Graphics::Draw3dBudding() {
  double m_rad = space_->radius;
  double d_rad = space_->bud_radius;
  double m_d_dist = space_->bud_height;
  glDisable(GL_LIGHTING);
  glDisable(GL_CULL_FACE);
  int lats = 12;
  int longs = 12;
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glLineWidth(1.0);
  double m_angle = acos((SQR(m_d_dist) + SQR(m_rad) - SQR(d_rad)) /
                        (2.0 * m_d_dist * m_rad));
  double d_angle = acos((SQR(m_d_dist) + SQR(d_rad) - SQR(m_rad)) /
                        (2.0 * m_d_dist * d_rad));
  for (int i = 0; i <= lats; i++) {
    double lat0 = -0.5 * M_PI + (M_PI - m_angle) * (double)(i - 1) / lats;
    double z0 = m_rad * sin(lat0) - z_correct_;
    double zr0 = cos(lat0);
    double lat1 = -0.5 * M_PI + (M_PI - m_angle) * (double)i / lats;
    double z1 = m_rad * sin(lat1) - z_correct_;
    double zr1 = cos(lat1);
    glBegin(GL_QUAD_STRIP);
    for (int j = 0; j <= longs; j++) {
      double lng = 2 * M_PI * (double)(j - 1) / longs;
      double x = m_rad * cos(lng);
      double y = m_rad * sin(lng);
      glNormal3f(x * zr0, y * zr0, z0);
      glVertex3f(x * zr0, y * zr0, z0);
      glNormal3f(x * zr1, y * zr1, z1);
      glVertex3f(x * zr1, y * zr1, z1);
    }
    glEnd();
  }
  lats = 10;
  longs = 10;
  for (int i = 0; i <= lats; i++) {
    double lat0 = d_angle - 0.5 * M_PI + (M_PI - d_angle) * (double)(i) / lats;
    double z0 = d_rad * sin(lat0) - z_correct_;
    double zr0 = cos(lat0);
    double lat1 =
        d_angle - 0.5 * M_PI + (M_PI - d_angle) * (double)(i + 1) / lats;
    double z1 = d_rad * sin(lat1) - z_correct_;
    double zr1 = cos(lat1);
    glBegin(GL_QUAD_STRIP);
    for (int j = 0; j <= longs; j++) {
      double lng = 2 * M_PI * (double)(j - 1) / longs;
      double x = d_rad * cos(lng);
      double y = d_rad * sin(lng);
      glNormal3f(x * zr0, y * zr0, z0 + m_d_dist);
      glVertex3f(x * zr0, y * zr0, z0 + m_d_dist);
      glNormal3f(x * zr1, y * zr1, z1 + m_d_dist);
      glVertex3f(x * zr1, y * zr1, z1 + m_d_dist);
    }
    glEnd();
  }

  // Turn off wireframe mode
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glEnable(GL_LIGHTING);
  glEnable(GL_CULL_FACE);
}

void Graphics::Draw2dBudding() {
  double m_rad = space_->radius;
  double d_rad = space_->bud_radius;
  double m_d_dist = space_->bud_height;
  // glUseProgram(0);

  /* Draw the bounding box. */
  // GLfloat box_color[4] = {1.0, 0.0, 0.5, 1.0};
  GLfloat box_color[4] = {0.5, 0.5, 0.5, 1.0};
  GLfloat v0, v1;
  GLfloat w = 4;
  glLineWidth(w);
  glColor4fv(box_color);
  glBegin(GL_LINE_LOOP);
  int n_vs = 100;
  double m_angle = acos((SQR(m_d_dist) + SQR(m_rad) - SQR(d_rad)) /
                        (2.0 * m_d_dist * m_rad));
  double d_angle = acos((SQR(m_d_dist) + SQR(d_rad) - SQR(m_rad)) /
                        (2.0 * m_d_dist * d_rad));
  double theta = 0.5 * M_PI + m_angle;
  double d_theta = (2.0 * M_PI - 2.0 * m_angle) / n_vs;
  double y_correct = 0.5 * (m_d_dist + m_rad + d_rad) - m_rad;
  for (int i = 0; i < n_vs; ++i) {
    v0 = m_rad * cos(theta + i * d_theta);
    v1 = m_rad * sin(theta + i * d_theta) - y_correct;
    glVertex2f(v0, v1);
  }
  theta = -0.5 * M_PI + d_angle;
  d_theta = (2.0 * M_PI - 2.0 * d_angle) / n_vs;
  for (int i = 0; i < n_vs; ++i) {
    v0 = d_rad * cos(theta + i * d_theta);
    v1 = m_d_dist + d_rad * sin(theta + i * d_theta) - y_correct;
    glVertex2f(v0, v1);
  }

  glEnd();
}

void Graphics::UpdateWindow() {
  /* Recalculate the view */
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  if (n_dim_ == 2) {
    glLoadIdentity();
    glTranslatef(xTrans_, yTrans_, 0.0);
  } else if (n_dim_ == 3) {
    glTranslatef(0.0, 0.0, zTrans_);
    glRotatef(xAngle_, 1.0, 0.0, 0.0);
    glRotatef(yAngle_, 0.0, 1.0, 0.0);
    glRotatef(zAngle_, 0.0, 0.0, 1.0);
  }
  glScalef(xyzScale_, xyzScale_, xyzScale_);

  /* Set square aspect ratio from current window size */
  // COMMENTED OUT TO GET RID OF UGLY 2D SNOWMAN - JMM
  // int width, height;
  // glfwGetFramebufferSize(window_, &width, &height);

  // int minsize = MIN(width, height);
  // glViewport((width - minsize)/2, (height - minsize)/2, minsize, minsize);

  /* Clear background */
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

/* Make a gl buffer */
GLuint GraphicsPrimitive::MakeBuffer(GLenum target, const void *buffer_data,
                                     GLsizei buffer_size) {
  GLuint buffer;
  glGenBuffers(1, &buffer);
  glBindBuffer(target, buffer);
  glBufferData(target, buffer_size, buffer_data, GL_STATIC_DRAW);
  return buffer;
}

void GraphicsPrimitive::MakeVertexBuffer() {
  vertex_buffer_ = MakeBuffer(GL_ARRAY_BUFFER, vertex_buffer_data_,
                              sizeof(GLfloat) * n_coords_);
}

void GraphicsPrimitive::MakeElementBuffer() {
  element_buffer_ = MakeBuffer(GL_ELEMENT_ARRAY_BUFFER, element_buffer_data_,
                               sizeof(GLushort) * n_triangles_ * 3);
}

/* Displays what part of our vertex shader failed to compile. This should never
   be used since the shader is hardcoded. */
void GraphicsPrimitive::ShowInfoLog(GLuint object,
                                    PFNGLGETSHADERIVPROC glGet__iv,
                                    PFNGLGETSHADERINFOLOGPROC glGet__InfoLog) {
  GLint log_length;
  char *log;

  glGet__iv(object, GL_INFO_LOG_LENGTH, &log_length);
  log = new char[log_length];
  glGet__InfoLog(object, log_length, NULL, log);
  fprintf(stderr, "%s", log);
  delete log;
}

void Graphics::InitDiscoRectangle() {
  discorectangle_.n_coords_ = (9 + 1) * 4;
  discorectangle_.n_triangles_ = (8 * 2 + 2);
  discorectangle_.vertex_buffer_data_ = new GLfloat[discorectangle_.n_coords_];
  discorectangle_.element_buffer_data_ =
      new GLushort[discorectangle_.n_triangles_ * 3];
  int i;
  GLfloat dphi = M_PI / 8.0;
  GLfloat phi;
  int index = 0;
  discorectangle_.vertex_buffer_data_[index + 0] = 0.0;
  discorectangle_.vertex_buffer_data_[index + 1] = 1e-6;
  index += 2;

  phi = 0.0;
  for (i = 0; i < 9; ++i) {
    discorectangle_.vertex_buffer_data_[index + 0] = 0.5 * cos(phi);
    discorectangle_.vertex_buffer_data_[index + 1] = 0.5 * sin(phi) + 1e-6;
    index += 2;
    phi += dphi;
  }

  int k = 0;
  for (i = 1; i < 9; ++i) {
    discorectangle_.element_buffer_data_[k] = 0;
    discorectangle_.element_buffer_data_[k + 1] = i;
    discorectangle_.element_buffer_data_[k + 2] = i + 1;
    k += 3;
  }

  discorectangle_.vertex_buffer_data_[index + 0] = 0.0;
  discorectangle_.vertex_buffer_data_[index + 1] = -1e-6;

  int center_vert = index / 2;
  index += 2;
  phi = M_PI;
  for (i = 0; i < 9; ++i) {
    discorectangle_.vertex_buffer_data_[index + 0] = 0.5 * cos(phi);
    discorectangle_.vertex_buffer_data_[index + 1] = 0.5 * sin(phi) - 1e-6;
    index += 2;
    phi += dphi;
  }

  for (i = 1; i < 9; ++i) {
    discorectangle_.element_buffer_data_[k] = center_vert;
    discorectangle_.element_buffer_data_[k + 1] = center_vert + i;
    discorectangle_.element_buffer_data_[k + 2] = center_vert + i + 1;
    k += 3;
  }

  discorectangle_.element_buffer_data_[k] = 1;
  discorectangle_.element_buffer_data_[k + 1] = center_vert + 1;
  discorectangle_.element_buffer_data_[k + 2] = center_vert + 9;
  k += 3;
  discorectangle_.element_buffer_data_[k] = 9;
  discorectangle_.element_buffer_data_[k + 1] = center_vert + 1;
  discorectangle_.element_buffer_data_[k + 2] = 1;
  k += 3;

  discorectangle_.MakeVertexBuffer();
  discorectangle_.MakeElementBuffer();

  /* Vertex shader source code */
  const char *vs =
      "#version 120\n"
      "attribute vec2 position;"
      "uniform float half_l;"
      "uniform float diameter;"
      "void main()"
      "{"
      "  vec4 stretched_pos;"
      "  if (position.y >= 0.0)"
      "    stretched_pos = vec4(diameter * position,0.0,1.0) + "
      "vec4(0.0,half_l,0.0,0.0);"
      "  else"
      "    stretched_pos = vec4(diameter * position,0.0,1.0) - "
      "vec4(0.0,half_l,0.0,0.0);"
      "  gl_FrontColor = gl_Color;"
      "  gl_Position = gl_ModelViewProjectionMatrix * (stretched_pos);"
      "} "
      "";

  GLint shader_ok;
  discorectangle_.vertex_shader_ = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(discorectangle_.vertex_shader_, 1, &vs, NULL);
  glCompileShader(discorectangle_.vertex_shader_);
  glGetShaderiv(discorectangle_.vertex_shader_, GL_COMPILE_STATUS, &shader_ok);
  if (!shader_ok) {
    fprintf(stderr, "Failed to compile vertex shader:\n");
    discorectangle_.ShowInfoLog(discorectangle_.vertex_shader_, glGetShaderiv,
                                glGetShaderInfoLog);
    glDeleteShader(discorectangle_.vertex_shader_);
    Logger::Error("Failed to compile vertex shader");
  }
  if (discorectangle_.vertex_shader_ == 0) {
    Logger::Error("Failed to compile vertex shader (vertex_shader_ = 0)");
  }

  /* Fragment shader source code */
  const char *fs = "#version 120\n"
                   "void main()"
                   "{"
                   "gl_FragColor = gl_Color;"
                   "}";

  discorectangle_.fragment_shader_ = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(discorectangle_.fragment_shader_, 1, &fs, NULL);
  glCompileShader(discorectangle_.fragment_shader_);
  glGetShaderiv(discorectangle_.fragment_shader_, GL_COMPILE_STATUS,
                &shader_ok);
  if (!shader_ok) {
    fprintf(stderr, "Failed to compile fragment shader:\n");
    discorectangle_.ShowInfoLog(discorectangle_.fragment_shader_, glGetShaderiv,
                                glGetShaderInfoLog);
    glDeleteShader(discorectangle_.fragment_shader_);
    Logger::Error("Failed to compile fragment shader");
  }
  if (discorectangle_.fragment_shader_ == 0) {
    Logger::Error("Failed to compile fragment shader (fragment_shader_ = 0)");
  }

  discorectangle_.MakeProgram();

  /* Get integer locations of relevant shader parameters */
  discorectangle_.uniforms_.half_l =
      glGetUniformLocation(discorectangle_.program_, "half_l");
  discorectangle_.uniforms_.diameter =
      glGetUniformLocation(discorectangle_.program_, "diameter");
  discorectangle_.attributes_.position =
      glGetAttribLocation(discorectangle_.program_, "position");
}

void Graphics::InitSpheroCylinder() {
  /* Create appropriate sized buffer objects. */
  spherocylinder_.n_coords_ = (16 * 16 + 2) * 3;
  spherocylinder_.n_triangles_ = (32 * 16);
  spherocylinder_.vertex_buffer_data_ = new GLfloat[spherocylinder_.n_coords_];
  spherocylinder_.element_buffer_data_ =
      new GLushort[spherocylinder_.n_triangles_ * 3];
  int i, j;
  double theta, phi;
  double dtheta = M_PI / 16.0;
  double dphi = 2.0 * M_PI / 16.0;
  double x, y, z;
  int index = 3;

  spherocylinder_.vertex_buffer_data_[0] = 0.0;
  spherocylinder_.vertex_buffer_data_[1] = 0.0;
  spherocylinder_.vertex_buffer_data_[2] = 0.5;

  theta = dtheta;
  phi = 0.0;

  for (i = 0; i < 8; ++i) {
    phi = 0.0;
    z = 0.5 * cos(theta) + 1e-6;
    for (j = 0; j < 16; ++j) {
      x = 0.5 * cos(phi) * sin(theta);
      y = 0.5 * sin(phi) * sin(theta);

      spherocylinder_.vertex_buffer_data_[index + 0] = x;
      spherocylinder_.vertex_buffer_data_[index + 1] = y;
      spherocylinder_.vertex_buffer_data_[index + 2] = z;
      phi += dphi;
      index += 3;
    }
    theta += dtheta;
  }

  theta -= dtheta;
  for (i = 0; i < 8; ++i) {
    phi = 0.0;
    z = 0.5 * cos(theta) - 1e-6;
    for (j = 0; j < 16; ++j) {
      x = 0.5 * cos(phi) * sin(theta);
      y = 0.5 * sin(phi) * sin(theta);

      spherocylinder_.vertex_buffer_data_[index + 0] = x;
      spherocylinder_.vertex_buffer_data_[index + 1] = y;
      spherocylinder_.vertex_buffer_data_[index + 2] = z;

      phi += dphi;
      index += 3;
    }
    theta += dtheta;
  }

  spherocylinder_.vertex_buffer_data_[index] = 0.0;
  spherocylinder_.vertex_buffer_data_[index + 1] = 0.0;
  spherocylinder_.vertex_buffer_data_[index + 2] = -0.5;
  index += 3;

  int offset = 0;
  int k = 0;
  for (i = 0; i < 15; ++i) {
    spherocylinder_.element_buffer_data_[k] = 0;
    spherocylinder_.element_buffer_data_[k + 1] = 1 + i;
    spherocylinder_.element_buffer_data_[k + 2] = 2 + i;
    k += 3;
  }
  spherocylinder_.element_buffer_data_[k] = 0;
  spherocylinder_.element_buffer_data_[k + 1] = 16;
  spherocylinder_.element_buffer_data_[k + 2] = 1;
  k += 3;

  for (j = 0; j < 15; ++j) {
    for (i = 1; i < 16; ++i) {
      spherocylinder_.element_buffer_data_[k] = i + j * 16;
      spherocylinder_.element_buffer_data_[k + 1] = i + 16 + j * 16;
      spherocylinder_.element_buffer_data_[k + 2] = i + 17 + j * 16;
      k += 3;
      spherocylinder_.element_buffer_data_[k] = i + 1 + j * 16;
      spherocylinder_.element_buffer_data_[k + 1] = i + j * 16;
      spherocylinder_.element_buffer_data_[k + 2] = i + 17 + j * 16;
      k += 3;
    }
    i = 16;
    spherocylinder_.element_buffer_data_[k] = i + j * 16;
    spherocylinder_.element_buffer_data_[k + 1] = i + 16 + j * 16;
    spherocylinder_.element_buffer_data_[k + 2] = i + 1 + j * 16;
    k += 3;
    spherocylinder_.element_buffer_data_[k] = i - 15 + j * 16;
    spherocylinder_.element_buffer_data_[k + 1] = i + j * 16;
    spherocylinder_.element_buffer_data_[k + 2] = i + 1 + j * 16;
    k += 3;
  }

  offset = spherocylinder_.n_coords_ / 3 - 17;
  for (i = offset; i < offset + 15; ++i) {
    spherocylinder_.element_buffer_data_[k] = spherocylinder_.n_coords_ / 3 - 1;
    spherocylinder_.element_buffer_data_[k + 1] = 1 + i;
    spherocylinder_.element_buffer_data_[k + 2] = 0 + i;

    k += 3;
  }
  spherocylinder_.element_buffer_data_[k] = spherocylinder_.n_coords_ / 3 - 1;
  spherocylinder_.element_buffer_data_[k + 1] = offset;
  spherocylinder_.element_buffer_data_[k + 2] =
      spherocylinder_.n_coords_ / 3 - 2;

  spherocylinder_.MakeVertexBuffer();

  /* It seems that just using the vertex buffer for normal information (since
    this starts out as a sphere) makes a much prettier spherocylinder */
  spherocylinder_.MakeElementBuffer();

  /* Vertex shader source code */
  const char *vs =
      "#version 120\n"
      "attribute vec4 position;"
      "uniform float half_l;"
      "uniform float diameter;"
      "varying vec4 diffuse,ambientGlobal,ambient;"
      "varying vec3 normal,lightDir,halfVector;"
      "varying float dist;"
      "void main()"
      "{"
      "  vec4 ecPos;"
      "  vec3 aux;"
      "  vec4 stretched_pos;"
      "  normal = normalize(gl_NormalMatrix * position.xyz);"
      "  if (position.z >= 0.0)"
      "    stretched_pos = diameter*position + vec4(0.0,0.0,half_l,0.0);"
      "  else"
      "    stretched_pos = diameter*position - vec4(0.0,0.0,half_l,0.0);"
      "  stretched_pos.w = position.w;"
      "  /* these are the new lines of code to compute the light's direction */"
      "  ecPos = gl_ModelViewMatrix * (stretched_pos);"
      "  aux = vec3(gl_LightSource[0].position-ecPos);"
      "  lightDir = normalize(aux);"
      "  dist = length(aux);"
      "  halfVector = normalize(gl_LightSource[0].halfVector.xyz);"
      ""
      "  /* Compute the diffuse, ambient and globalAmbient terms */"
      "  diffuse = gl_FrontMaterial.diffuse * gl_LightSource[0].diffuse;"
      ""
      "  /* The ambient terms have been separated since one of them */"
      "  /* suffers attenuation */"
      "  ambient = gl_FrontMaterial.ambient * gl_LightSource[0].ambient;"
      "  ambientGlobal = gl_LightModel.ambient * gl_FrontMaterial.ambient;"
      ""
      "  gl_Position = gl_ModelViewProjectionMatrix * (stretched_pos);"
      "} "
      "";

  /* Build vertex shader from source specified in vs string */
  GLint shader_ok;
  spherocylinder_.vertex_shader_ = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(spherocylinder_.vertex_shader_, 1, &vs, NULL);
  glCompileShader(spherocylinder_.vertex_shader_);
  glGetShaderiv(spherocylinder_.vertex_shader_, GL_COMPILE_STATUS, &shader_ok);
  if (!shader_ok) {
    fprintf(stderr, "Failed to compile vertex shader:\n");
    spherocylinder_.ShowInfoLog(spherocylinder_.vertex_shader_, glGetShaderiv,
                                glGetShaderInfoLog);
    glDeleteShader(spherocylinder_.vertex_shader_);
    Logger::Error("Failed to compile vertex shader");
  }
  if (spherocylinder_.vertex_shader_ == 0) {
    Logger::Error("Failed to compile vertex shader (vertex_shader_ = 0)");
  }
  /* Fragment shader source code */
  std::ostringstream ss;
  ss << "#version 120\n"
        "varying vec4 diffuse,ambientGlobal, ambient;"
        "varying vec3 normal,lightDir,halfVector;"
        "varying float dist;"
        "void main()"
        "{"
        "vec3 n,halfV,viewV,ldir;"
        "float NdotL,NdotHV;"
        "vec4 color = ambientGlobal;"
        "float att;"
        "/* a fragment shader can't write a varying variable, hence we need"
        "a new variable to store the normalized interpolated normal */"
        "n = normalize(normal);"
        ""
        "    /* compute the dot product between normal and normalized lightdir "
        "*/"
        "NdotL = max(dot(n,normalize(lightDir)),0.0);"
        ""
        "if (NdotL > 0.0) {"
        ""
        "att = 1.0 / (gl_LightSource[0].constantAttenuation +"
        "gl_LightSource[0].linearAttenuation * dist +"
        "gl_LightSource[0].quadraticAttenuation * dist * dist);"
        "color += att * (diffuse * NdotL + ambient);"
        "}"
        ""
        "color[3] = " //Set object opacity here
     << alpha_
     << ";"
        "gl_FragColor = color;"
        "}";

  std::string s = ss.str(); // Need an intermediary string for this to work
  const char *fs = s.c_str();

  /* Build Fragment Shader from source specified in fs string */
  spherocylinder_.fragment_shader_ = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(spherocylinder_.fragment_shader_, 1, &fs, NULL);
  glCompileShader(spherocylinder_.fragment_shader_);
  glGetShaderiv(spherocylinder_.fragment_shader_, GL_COMPILE_STATUS,
                &shader_ok);
  if (!shader_ok) {
    fprintf(stderr, "Failed to compile fragment shader:\n");
    spherocylinder_.ShowInfoLog(spherocylinder_.fragment_shader_, glGetShaderiv,
                                glGetShaderInfoLog);
    glDeleteShader(spherocylinder_.fragment_shader_);
    Logger::Error("Failed to compile fragment shader");
  }

  spherocylinder_.MakeProgram();

  /* Get integer locations of relevant shader parameters */
  spherocylinder_.uniforms_.half_l =
      glGetUniformLocation(spherocylinder_.program_, "half_l");
  spherocylinder_.uniforms_.diameter =
      glGetUniformLocation(spherocylinder_.program_, "diameter");
  spherocylinder_.attributes_.position =
      glGetAttribLocation(spherocylinder_.program_, "position");
}

/* Initialize our custom shader program for spherocylinders */
void GraphicsPrimitive::MakeProgram() {
  program_ = glCreateProgram();
  glAttachShader(program_, vertex_shader_);
  glAttachShader(program_, fragment_shader_);
  glLinkProgram(program_);

  GLint program_ok;
  glGetProgramiv(program_, GL_LINK_STATUS, &program_ok);
  if (!program_ok) {
    fprintf(stderr, "Failed to link shader program:\n");
    ShowInfoLog(program_, glGetProgramiv, glGetProgramInfoLog);
    glDeleteProgram(program_);
    Logger::Error("Failed to link shader program");
  }
}

// void GraphicsText::MakeProgram() {
// program_ = glCreateProgram();
////vertex_shader = create_shader(vertexfile, GL_VERTEX_SHADER);
// glAttachShader(program_, vertex_shader_);
////fragment_shader = create_shader(fragmentfile, GL_FRAGMENT_SHADER);
// glAttachShader(program_, fragment_shader_);
// glLinkProgram(program_);
// GLint link_ok = GL_FALSE;
// glGetProgramiv(program_, GL_LINK_STATUS, &link_ok);
// if (!link_ok) {
// fprintf(stderr, "Failed to link text shader program:\n");
// ShowInfoLog(program_, glGetProgramiv, glGetProgramInfoLog);
// glDeleteProgram(program_);
// exit(1);
//}
// if (program_ == 0) {
// fprintf(stderr, "Failed to create text program\n");
//}
//}

// void Graphics::InitText() {
// if (FT_Init_FreeType(&(text_.ft_))) {
// fprintf(stderr, "Could not init freetype library\n");
// exit(1);
//}
//[> Load a font <]
// if (FT_New_Face(text_.ft_,"FreeSans.ttf", 0, &(text_.face_))) {
// fprintf(stderr, "Could not open font FreeSans.ttf\n");
// exit(1);
//}

//[> Vertex shader source code <]
// const char *vs =
//"#version 120\n"
//"attribute vec4 coord;"
//"varying vec2 texpos;"
//"void main(void) {"
//"gl_Position = vec4(coord.xy, 0, 1);"
//"texpos = coord.zw;"
//"}";
//[> Build vertex shader from source specified in vs string <]
// GLint shader_ok;
// text_.vertex_shader_ = glCreateShader(GL_VERTEX_SHADER);
// glShaderSource(text_.vertex_shader_, 1, &vs, NULL);
// glCompileShader(text_.vertex_shader_);
// glGetShaderiv(text_.vertex_shader_, GL_COMPILE_STATUS, &shader_ok);
// if (!shader_ok) {
// fprintf(stderr, "Failed to compile text vertex shader:\n");
// text_.ShowInfoLog(text_.vertex_shader_, glGetShaderiv, glGetShaderInfoLog);
// glDeleteShader(text_.vertex_shader_);
// exit(1);
//}
// if (text_.vertex_shader_ == 0)
// exit(1);

// const char *fs =
//"#version 120\n"
//"varying vec2 texpos;"
//"uniform sampler2D tex;"
//"uniform vec4 color;"
//"void main() {"
//"gl_FragColor = vec4(1, 1, 1, texture2D(tex, texpos).a) * color;"
//"}";

//[> Build Fragment Shader from source specified in fs string <]
// text_.fragment_shader_ = glCreateShader(GL_FRAGMENT_SHADER);
// glShaderSource(text_.fragment_shader_, 1, &fs, NULL);
// glCompileShader(text_.fragment_shader_);
// glGetShaderiv(text_.fragment_shader_, GL_COMPILE_STATUS, &shader_ok);
// if (!shader_ok) {
// fprintf(stderr, "Failed to compile text fragment shader:\n");
// text_.ShowInfoLog(text_.fragment_shader_, glGetShaderiv, glGetShaderInfoLog);
// glDeleteShader(text_.fragment_shader_);
// exit(1);
//}

// text_.MakeProgram();

// text_.attribute_coord_ = glGetAttribLocation(text_.program_, "coord");
// text_.uniform_tex_ = glGetUniformLocation(text_.program_, "tex");
// text_.uniform_color_ = glGetUniformLocation(text_.program_, "color");
// if(text_.attribute_coord_ == -1 || text_.uniform_tex_ == -1 ||
// text_.uniform_color_ == -1) {
// fprintf(stderr, "Could not bind text uniform {%d %d
// %d}\n",text_.attribute_coord_, text_.uniform_tex_, text_.uniform_color_);
// exit(1);
//}
//// Create the vertex buffer object
// glGenBuffers(1, &(text_.buffer_));

//}

// void GraphicsText::ShowInfoLog(GLuint object,
// PFNGLGETSHADERIVPROC glGet__iv,
// PFNGLGETSHADERINFOLOGPROC glGet__InfoLog) {
// GLint log_length;
// char *log;

// glGet__iv(object, GL_INFO_LOG_LENGTH, &log_length);
// log = new char[log_length];
// glGet__InfoLog(object, log_length, NULL, log);
// fprintf(stderr, "%s", log);
// delete log;
//}

// void GraphicsText::RenderText(const char *text, float x, float y, int sx, int
// sy) { const char *p; FT_GlyphSlot g = face_->glyph;

//[> Create a texture that will be used to hold one "glyph" <]
// GLuint tex;

// glActiveTexture(GL_TEXTURE0);
// glGenTextures(1, &tex);
// glBindTexture(GL_TEXTURE_2D, tex);
// glUniform1i(uniform_tex_, 0);

//[> We require 1 byte alignment when uploading texture data <]
// glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

//[> Clamping to edges is important to prevent artifacts when scaling <]
// glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
// glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

//[> Linear filtering usually looks best for text <]
// glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
// glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

//[> Set up the VBO for our vertex data <]
// glEnableVertexAttribArray(attribute_coord_);
// glBindBuffer(GL_ARRAY_BUFFER, buffer_);
// glVertexAttribPointer(attribute_coord_, 4, GL_FLOAT, GL_FALSE, 0, 0);

//[> Loop through all characters <]
// for (p = text; *p; p++) {
//[> Try to load and render the character <]
// if (FT_Load_Char(face_, *p, FT_LOAD_RENDER))
// continue;

//[> Upload the "bitmap", which contains an 8-bit grayscale image, as an alpha
// texture <] glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA, g->bitmap.width,
// g->bitmap.rows, 0, GL_ALPHA, GL_UNSIGNED_BYTE, g->bitmap.buffer);

//[> Calculate the vertex and texture coordinates <]
// float x2 = x + g->bitmap_left * sx;
// float y2 = -y - g->bitmap_top * sy;
// float w = g->bitmap.width * sx;
// float h = g->bitmap.rows * sy;

// point box[4] = {
//{x2, -y2, 0, 0},
//{x2 + w, -y2, 1, 0},
//{x2, -y2 - h, 0, 1},
//{x2 + w, -y2 - h, 1, 1},
//};

//[> Draw the character on the screen <]
// glBufferData(GL_ARRAY_BUFFER, sizeof box, box, GL_DYNAMIC_DRAW);
// glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

//[> Advance the cursor to the start of the next character <]
// x += (g->advance.x >> 6) * sx;
// y += (g->advance.y >> 6) * sy;
//}

// glDisableVertexAttribArray(attribute_coord_);
// glDeleteTextures(1, &tex);
//}

#endif
