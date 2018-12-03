// g2o - General Graph Optimization
// Copyright (C) 2011 R. Kuemmerle, G. Grisetti, H. Strasdat, W. Burgard
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "edge_se3.h"
#include "isometry3d_gradients.h"
#include <iostream>

#ifdef G2O_HAVE_OPENGL
#include "g2o/stuff/opengl_wrapper.h"
#include "g2o/stuff/opengl_primitives.h"
#endif

namespace g2o {
  using namespace std;

  EdgeSE3::EdgeSE3() : BaseBinaryEdge<6, Isometry3, VertexSE3, VertexSE3>() {
    information().setIdentity();
  }

  bool EdgeSE3::read(std::istream& is) {
    Vector7 meas;
    for (int i=0; i<7; i++) 
      is >> meas[i];
    // normalize the quaternion to recover numerical precision lost by storing as human readable text
    Vector4::MapType(meas.data()+3).normalize();
    setMeasurement(internal::fromVectorQT(meas));

    if (is.bad()) {
      return false;
    }
    for ( int i=0; i<information().rows() && is.good(); i++)
      for (int j=i; j<information().cols() && is.good(); j++){
        is >> information()(i,j);
        if (i!=j)
          information()(j,i)=information()(i,j);
      }
    if (is.bad()) {
      //  we overwrite the information matrix with the Identity
      information().setIdentity();
    } 
    return true;
  }

  bool EdgeSE3::write(std::ostream& os) const {
    Vector7 meas=internal::toVectorQT(_measurement);
    for (int i=0; i<7; i++) os  << meas[i] << " ";
    for (int i=0; i<information().rows(); i++)
      for (int j=i; j<information().cols(); j++) {
        os <<  information()(i,j) << " ";
      }
    return os.good();
  }

  void EdgeSE3::computeError() {
    VertexSE3 *from = static_cast<VertexSE3*>(_vertices[0]);
    VertexSE3 *to   = static_cast<VertexSE3*>(_vertices[1]);
    Isometry3 delta=_inverseMeasurement * from->estimate().inverse() * to->estimate();
    _error=internal::toVectorMQT(delta);
  }

  bool EdgeSE3::setMeasurementFromState(){
    VertexSE3 *from = static_cast<VertexSE3*>(_vertices[0]);
    VertexSE3 *to   = static_cast<VertexSE3*>(_vertices[1]);
    Isometry3 delta = from->estimate().inverse() * to->estimate();
    setMeasurement(delta);
    return true;
  }
  
  void EdgeSE3::linearizeOplus(){
    
    // BaseBinaryEdge<6, Isometry3, VertexSE3, VertexSE3>::linearizeOplus();
    // return;

    VertexSE3 *from = static_cast<VertexSE3*>(_vertices[0]);
    VertexSE3 *to   = static_cast<VertexSE3*>(_vertices[1]);
    Isometry3 E;
    const Isometry3& Xi=from->estimate();
    const Isometry3& Xj=to->estimate();
    const Isometry3& Z=_measurement;
    internal::computeEdgeSE3Gradient(E, _jacobianOplusXi , _jacobianOplusXj, Z, Xi, Xj);
  }

  void EdgeSE3::initialEstimate(const OptimizableGraph::VertexSet& from_, OptimizableGraph::Vertex* /*to_*/) {
    VertexSE3 *from = static_cast<VertexSE3*>(_vertices[0]);
    VertexSE3 *to   = static_cast<VertexSE3*>(_vertices[1]);

    if (from_.count(from) > 0) {
      to->setEstimate(from->estimate() * _measurement);
    } else
      from->setEstimate(to->estimate() * _measurement.inverse());
    //cerr << "IE" << endl;
  }

  EdgeSE3WriteGnuplotAction::EdgeSE3WriteGnuplotAction(): WriteGnuplotAction(typeid(EdgeSE3).name()){}

  HyperGraphElementAction* EdgeSE3WriteGnuplotAction::operator()(HyperGraph::HyperGraphElement* element, HyperGraphElementAction::Parameters* params_){
    if (typeid(*element).name()!=_typeName)
      return nullptr;
    WriteGnuplotAction::Parameters* params=static_cast<WriteGnuplotAction::Parameters*>(params_);
    if (!params->os){
      std::cerr << __PRETTY_FUNCTION__ << ": warning, on valid os specified" << std::endl;
      return nullptr;
    }

    EdgeSE3* e =  static_cast<EdgeSE3*>(element);
    VertexSE3* fromEdge = static_cast<VertexSE3*>(e->vertices()[0]);
    VertexSE3* toEdge   = static_cast<VertexSE3*>(e->vertices()[1]);
    Vector6 fromV, toV;
    fromV=internal::toVectorMQT(fromEdge->estimate());
    toV=internal::toVectorMQT(toEdge->estimate());
    for (int i=0; i<6; i++){
      *(params->os) << fromV[i] << " ";
    }
    for (int i=0; i<6; i++){
      *(params->os) << toV[i] << " ";
    }
    *(params->os) << std::endl;
    return this;
  }

#ifdef G2O_HAVE_OPENGL
  bool EdgeSE3DrawAction::refreshPropertyPtrs(HyperGraphElementAction::Parameters* params_){
    if (!DrawAction::refreshPropertyPtrs(params_))
      return false;
    if (_previousParams){
      _showMeasurementAndError = _previousParams->makeProperty<BoolProperty>(_typeName + "::SHOW_MEASUREMENT_AND_ERROR", false);
      _showEllipsoid = _previousParams->makeProperty<BoolProperty>(_typeName + "::SHOW_STD_DEV", false);
    } else {
      _showMeasurementAndError = 0;
      _showEllipsoid = 0;
    }
    return true;
  }

  EdgeSE3DrawAction::EdgeSE3DrawAction(): DrawAction(typeid(EdgeSE3).name()){}

  void  EdgeSE3DrawAction::drawMeasurementAndError(Eigen::Vector3f& fromPos,
                                                   Eigen::Vector3f& estToPos,
                                                   Eigen::Vector3f& measToPos)
  {
      glPushAttrib(GL_LINE_BIT);
      glLineStipple(1, 0xAAAA);
      glLineWidth(EDGE_LINE_WIDTH);
      glEnable(GL_LINE_STIPPLE);
      glBegin(GL_LINES);
      //Measured transformation in yellow
      glColor3f(POSE_EDGE_MEASUREMENT_COLOR);
      glVertex3f(fromPos.x(),fromPos.y(),fromPos.z());
      glVertex3f(measToPos.x(),measToPos.y(),measToPos.z());
      //and difference to estimate in dotted red
      glColor3f(POSE_EDGE_ERROR_COLOR);
      glVertex3f(measToPos.x(),measToPos.y(),measToPos.z());
      glVertex3f(estToPos.x(),estToPos.y(),estToPos.z());
      glEnd();
      glPopAttrib();
  }

  void EdgeSE3DrawAction::drawUncertainty(Isometry3& measuredTo, EdgeSE3::InformationType& infoMat)
  {
      //Draw uncertainty ellipsoid for one std dev
      glColor3f(EDGE_UNCERTAINTY_ELLIPSOID_COLOR);
      glPushMatrix();
      glPushAttrib(GL_POLYGON_BIT);
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);//Draw wireframe
      glMultMatrixd(measuredTo.matrix().data());
      EdgeSE3::InformationType cov = infoMat.inverse();
      opengl::drawEllipsoid(3*sqrt(cov(0,0)), 3*sqrt(cov(1,1)), 3*sqrt(cov(2,2)));
      glPopAttrib();//Restore from wireframe
      glPopMatrix();
  }

  HyperGraphElementAction* EdgeSE3DrawAction::operator()(HyperGraph::HyperGraphElement* element, 
               HyperGraphElementAction::Parameters* params_){
    if (typeid(*element).name()!=_typeName)
      return nullptr;
    refreshPropertyPtrs(params_);
    if (! _previousParams)
      return this;
    
    if (_show && !_show->value())
      return this;
    
    auto* edge = static_cast<EdgeSE3*>(element);
    auto* from = static_cast<VertexSE3*>(edge->vertices()[0]);
    auto* to   = static_cast<VertexSE3*>(edge->vertices()[1]);
    Eigen::Vector3f fromPos = from->estimate().translation().cast<float>();
    Eigen::Vector3f estToPos = to->estimate().translation().cast<float>();
    Isometry3 measuredTo = (from->estimate() * edge->measurement());
    Eigen::Vector3f measToPos = measuredTo.translation().cast<float>();
    if (! from || ! to)
      return this;
    glPushAttrib(GL_ENABLE_BIT);
    glPushAttrib(GL_LINE_BIT);
    glDisable(GL_LIGHTING);
    glLineWidth(EDGE_LINE_WIDTH);
    glBegin(GL_LINES);
    glColor3f(POSE_EDGE_COLOR);
    glVertex3f(fromPos.x(),fromPos.y(),fromPos.z());
    glVertex3f(estToPos.x(),estToPos.y(),estToPos.z());
    glEnd();

    if(_showMeasurementAndError && _showMeasurementAndError->value()){
      drawMeasurementAndError(fromPos, estToPos, measToPos);
    }
     if(_showEllipsoid && _showEllipsoid->value()){
      drawUncertainty(measuredTo, edge->information());
    }

    glPopAttrib();//restore Line width
    glPopAttrib();//restore enable bit (lighting?)
    return this;
  }
#endif

}
