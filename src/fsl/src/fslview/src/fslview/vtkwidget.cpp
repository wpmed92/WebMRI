/*  FSLView - 2D/3D Interactive Image Viewer

    Authors:    Rama Aravind Vorray
		James Saunders
		David Flitney 
		Mark Jenkinson
		Stephen Smith

    FMRIB Image Analysis Group

    Copyright (C) 2002-2003 University of Oxford  */

/*  CCOPYRIGHT */

#if defined(WIN32)
#pragma warning(disable:4786)
#endif

#include <vector>
#include <map>

#include <qapplication.h>
#include <qfiledialog.h>

#include "vtkwidget.h"
#include "maintoolbar.h"
#include "metaimage.h"
#include "tracker.h"

#include <vtkImageActor.h>
#include <vtkVolume.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkLookupTable.h>
#include <vtkImageBlend.h>
#include <vtkImagePlaneWidget.h>
#include <vtkImageMapToColors.h>
#include <vtkImageOpenClose3D.h>
#include <vtkVolumeRayCastCompositeFunction.h>
#include <vtkVolumeRayCastMapper.h>
#include <vtkColorTransferFunction.h>
#include <vtkMarchingCubes.h>
#include <vtkDiscreteMarchingCubes.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkClipPolyData.h>
#include <vtkClipVolume.h>
#include <vtkCutter.h>
#include <vtkDecimatePro.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataWriter.h>
#include <vtkDataSetWriter.h>
#include <vtkBoxWidget.h>
#include <vtkPlaneSource.h>
#include <vtkPlanes.h>
#include <vtkPlane.h>
#include <vtkImageThreshold.h>
#include <vtkImageReslice.h>
#include <vtkProperty.h>
#include <vtkVolumeProperty.h>
#include <vtkPiecewiseFunction.h>
#include <vtkScalarsToColors.h>
#include <vtkTextureMapToPlane.h>
#include <vtkTransform.h>
#include <vtkDataSetMapper.h>
#include <vtkImageMapper.h>
#include <vtkProbeFilter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkImageStencil.h>
#include <vtkImplicitDataSet.h>
#include <vtkDataSetWriter.h>
#include <QVTKWidget.h>

#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkDataSetAttributes.h>
#include <vtkPolyDataNormals.h>
#include <vtkTensorGlyph.h>
#include <vtkSphereSource.h>

#include <vtkImageGaussianSmooth.h>

#include <vtkSmoothPolyDataFilter.h>
#include <vtkImageDilateErode3D.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkDepthSortPolyData.h>

#include <vtkWindowToImageFilter.h>
#include <vtkTIFFWriter.h>

#include "vtktoolbar.h"
#include "vtkpropertydialog.h"

class VTKLookUpTable
{
public:
  static vtkLookupTable* RedYellow() 
  {
    vtkPiecewiseFunction *r = vtkPiecewiseFunction::New();
    vtkPiecewiseFunction *g = vtkPiecewiseFunction::New();
    vtkPiecewiseFunction *b = vtkPiecewiseFunction::New();
    vtkPiecewiseFunction *a = vtkPiecewiseFunction::New();

    r->AddPoint(500,   1);
    r->AddPoint(32000, 1);
    g->AddPoint(500,   0);
    g->AddPoint(32000, 1);
    b->AddPoint(500,   0);
    b->AddPoint(32000, 0);
    a->AddPoint(500,   1);
    a->AddPoint(32000, 1);
    r->ClampingOff();
    g->ClampingOff();
    b->ClampingOff();
    a->ClampingOff();

    vtkLookupTable* ctfun = vtkLookupTable::New();
    ctfun->SetNumberOfTableValues(32001);    
    for(unsigned short i = 0; i < 32000; ++i)
      ctfun->SetTableValue(int(i),
			   r->GetValue(i), g->GetValue(i), b->GetValue(i), 
			   a->GetValue(i));
 
    ctfun->SetTableRange(0, 32000);
    ctfun->Build();

    return ctfun;
  }

  static vtkLookupTable* BlueLightblue()
  {
    vtkPiecewiseFunction *r = vtkPiecewiseFunction::New();
    vtkPiecewiseFunction *g = vtkPiecewiseFunction::New();
    vtkPiecewiseFunction *b = vtkPiecewiseFunction::New();
    vtkPiecewiseFunction *a = vtkPiecewiseFunction::New();

    r->AddPoint(500,   0);
    r->AddPoint(32000, 0);
    g->AddPoint(500,   0);
    g->AddPoint(32000, 1);
    b->AddPoint(500,   1);
    b->AddPoint(32000, 1);
    a->AddPoint(500,   1);
    a->AddPoint(32000, 1);
    r->ClampingOff();
    g->ClampingOff();
    b->ClampingOff();
    a->ClampingOff();

    vtkLookupTable* ctfun = vtkLookupTable::New();
    ctfun->SetNumberOfTableValues(32001);    
    for(unsigned short i = 0; i < 32000; ++i)
      ctfun->SetTableValue(int(i),
			   r->GetValue(i), g->GetValue(i), b->GetValue(i), 
			   a->GetValue(i));

    ctfun->SetTableRange(0, 32000);
    ctfun->Build();

    return ctfun;
  }

  static vtkLookupTable* GreyScale()
  {
    vtkPiecewiseFunction *r = vtkPiecewiseFunction::New();
    vtkPiecewiseFunction *g = vtkPiecewiseFunction::New();
    vtkPiecewiseFunction *b = vtkPiecewiseFunction::New();
    vtkPiecewiseFunction *a = vtkPiecewiseFunction::New();

    r->AddPoint(500,   0);
    r->AddPoint(32000, 1);
    g->AddPoint(500,   0);
    g->AddPoint(32000, 1);
    b->AddPoint(500,   0);
    b->AddPoint(32000, 1);
    a->AddPoint(500,   1);
    a->AddPoint(32000, 1);
    r->ClampingOff();
    g->ClampingOff();
    b->ClampingOff();
    a->ClampingOff();

    vtkLookupTable* ctfun = vtkLookupTable::New();
    ctfun->SetNumberOfTableValues(32001);    
    for(unsigned short i = 0; i < 32000; ++i)
      ctfun->SetTableValue(int(i),
			   r->GetValue(i), g->GetValue(i), b->GetValue(i), 
			   a->GetValue(i));

    ctfun->SetTableRange(0, 32000);
    ctfun->Build();
    
    return ctfun;
  }

  static vtkLookupTable* Red()
  {
    vtkPiecewiseFunction *r = vtkPiecewiseFunction::New();
    vtkPiecewiseFunction *g = vtkPiecewiseFunction::New();
    vtkPiecewiseFunction *b = vtkPiecewiseFunction::New();
    vtkPiecewiseFunction *a = vtkPiecewiseFunction::New();

    r->AddPoint(500,   1);
    r->AddPoint(32000, 1);
    g->AddPoint(500,   0);
    g->AddPoint(32000, 0);
    b->AddPoint(500,   0);
    b->AddPoint(32000, 0);
    a->AddPoint(500,   1);
    a->AddPoint(32000, 1);
    r->ClampingOff();
    g->ClampingOff();
    b->ClampingOff();
    a->ClampingOff();

    vtkLookupTable* ctfun = vtkLookupTable::New();
    ctfun->SetNumberOfTableValues(32001);
    for(unsigned short i = 0; i < 32000; ++i)
      ctfun->SetTableValue(int(i),
			   r->GetValue(i), g->GetValue(i), b->GetValue(i), 
			   a->GetValue(i));

    ctfun->SetTableRange(0, 32000);
    ctfun->Build();

    return ctfun;
  }

  static vtkLookupTable* Green()
  {
    vtkPiecewiseFunction *r = vtkPiecewiseFunction::New();
    vtkPiecewiseFunction *g = vtkPiecewiseFunction::New();
    vtkPiecewiseFunction *b = vtkPiecewiseFunction::New();
    vtkPiecewiseFunction *a = vtkPiecewiseFunction::New();

    r->AddPoint(500,   0);
    r->AddPoint(32000, 0);
    g->AddPoint(500,   1);
    g->AddPoint(32000, 1);
    b->AddPoint(500,   0);
    b->AddPoint(32000, 0);
    a->AddPoint(500,   1);
    a->AddPoint(32000, 1);
    r->ClampingOff();
    g->ClampingOff();
    b->ClampingOff();
    a->ClampingOff();

    vtkLookupTable* ctfun = vtkLookupTable::New();
    ctfun->SetNumberOfTableValues(32001);
    for(unsigned short i = 0; i < 32000; ++i)
      ctfun->SetTableValue(int(i),
			   r->GetValue(i), g->GetValue(i), b->GetValue(i), 
			   a->GetValue(i));

    ctfun->SetTableRange(0, 32000);
    ctfun->Build();

    return ctfun;
  }
  
  static vtkLookupTable* Blue()
  {
    vtkPiecewiseFunction *r = vtkPiecewiseFunction::New();
    vtkPiecewiseFunction *g = vtkPiecewiseFunction::New();
    vtkPiecewiseFunction *b = vtkPiecewiseFunction::New();
    vtkPiecewiseFunction *a = vtkPiecewiseFunction::New();

    r->AddPoint(500,   0);
    r->AddPoint(32000, 0);
    g->AddPoint(500,   0);
    g->AddPoint(32000, 0);
    b->AddPoint(500,   1);
    b->AddPoint(32000, 1);
    a->AddPoint(500,   1);
    a->AddPoint(32000, 1);
    r->ClampingOff();
    g->ClampingOff();
    b->ClampingOff();
    a->ClampingOff();

    vtkLookupTable* ctfun = vtkLookupTable::New();
    ctfun->SetNumberOfTableValues(32001);
    for(unsigned short i = 0; i < 32000; ++i)
      ctfun->SetTableValue(int(i),
			   r->GetValue(i), g->GetValue(i), b->GetValue(i), 
			   a->GetValue(i));

    ctfun->SetTableRange(0, 32000);
    ctfun->Build();

    return ctfun;
  }

  vtkLookupTable* GetNextLookupTable() const
  {
    m_index = m_index % m_luts.size();
    return m_luts.at(m_index++);
  }

  vtkLookupTable* GetLookupTableByName(const std::string& name)
  {
    if(name == "Greyscale")      return GreyScale();
    if(name == "Red-Yellow")     return RedYellow();
    if(name == "Blue-Lightblue") return BlueLightblue();
    if(name == "Red")            return Red();
    if(name == "Green")          return Green();
    if(name == "Blue")           return Blue();
    return GreyScale();
  }

  VTKLookUpTable(): m_index(1) 
  {
    m_luts.reserve(6);
    m_luts.push_back(GreyScale());
    m_luts.push_back(RedYellow());
    m_luts.push_back(BlueLightblue());
    m_luts.push_back(Red());
    m_luts.push_back(Blue());
    m_luts.push_back(Green());
  }

private:
  mutable unsigned int m_index;
  
  std::vector<vtkLookupTable*> m_luts;
};

class InteractorCallback: public vtkCommand
{
public:
  static InteractorCallback* New()
  { return new InteractorCallback; }
  virtual void Execute(vtkObject* caller, unsigned long eventid, void*)
  {
    cout << "An event!" << endl;
    switch(eventid) {
    case vtkCommand::StartInteractionEvent:
      cout << "Start event" << endl;
      if(m_props->inqOpacity() != 1.0)
	{
	  cout << "Disconnect sorter" << endl;
	  m_mapper->SetInput(m_clipper->GetOutput());
	}
      break;
    case vtkCommand::EndInteractionEvent:
      cout << "End event" << endl;
      if(m_props->inqOpacity() != 1.0)
	{
	  cout << "Connect sorter" << endl;
	  m_mapper->SetInput(m_sorter->GetOutput());
	}
      break;
    default:
      break;
    }
  }
  void SetProperties(VTKProperties* p) { m_props = p; }
  void SetClipper(vtkClipPolyData* c) { m_clipper = c; }
  void SetSorter(vtkDepthSortPolyData* s) { m_sorter = s; }
  void SetMapper(vtkPolyDataMapper* m) { m_mapper = m; }

private:
  
  vtkDepthSortPolyData* m_sorter;
  vtkClipPolyData* m_clipper;
  vtkPolyDataMapper* m_mapper;
  VTKProperties* m_props;
};

class ClippingBoxCallback: public vtkCommand
{
public:
  static ClippingBoxCallback* New()
  { return new ClippingBoxCallback; }
  virtual void Execute(vtkObject* caller, unsigned long eventid, void*)
  {
    vtkBoxWidget* box = reinterpret_cast<vtkBoxWidget*>(caller);
    box->GetPlanes(m_planes);
    box->GetTransform(m_xform);
  }
  void SetTransform(vtkTransform* t) { m_xform = t; }
  void SetPlanes(vtkPlanes* p) { m_planes = p; }

private:
  
  vtkPlanes* m_planes;
  vtkTransform* m_xform;
};

VTKProperties::VTKProperties(): 
  m_lowerThreshold(10), m_upperThreshold(32000), 
  m_mcThreshold(0.001),
  m_interpMode(VTK_GOURAUD),
  m_kernelX(1), m_kernelY(1), m_kernelZ(1),
  m_relaxationFactor(0.8), 
  m_iterations(5),
  m_ambient(0.1), m_diffuse(0.9), m_opacity(1.0), m_specular(0.9), m_specularPower(35),
  m_featureAngle(140.0),
  m_colorR(1.0), m_colorG(0.80), m_colorB(0.70)
{ 
  m_stdDev[0]=0.5;
  m_stdDev[1]=0.5;
  m_stdDev[2]=0.5;
  m_radii[0]=1.0;
  m_radii[1]=1.0;
  m_radii[2]=1.0;
} 

VTKProperties::VTKProperties(const VTKProperties& rhs): 
  m_lowerThreshold(rhs.m_lowerThreshold),
  m_upperThreshold(rhs.m_upperThreshold),
  m_mcThreshold(rhs.m_mcThreshold),
  m_interpMode(rhs.m_interpMode),
  m_kernelX(rhs.m_kernelX),
  m_kernelY(rhs.m_kernelY),
  m_kernelZ(rhs.m_kernelZ),
  m_relaxationFactor(rhs.m_relaxationFactor),
  m_iterations(rhs.m_iterations),
  m_ambient(rhs.m_ambient),
  m_diffuse(rhs.m_diffuse),
  m_opacity(rhs.m_opacity),
  m_specular(rhs.m_specular),
  m_specularPower(rhs.m_specularPower),
  m_featureAngle(rhs.m_featureAngle),
  m_colorR(rhs.m_colorR), m_colorG(rhs.m_colorG), m_colorB(rhs.m_colorB)
{
  m_stdDev[0]=rhs.m_stdDev[0];
  m_stdDev[1]=rhs.m_stdDev[1];
  m_stdDev[2]=rhs.m_stdDev[2];
  m_radii[0]=rhs.m_radii[0];
  m_radii[1]=rhs.m_radii[1];
  m_radii[2]=rhs.m_radii[2];
}

VTKProperties& VTKProperties::operator=(const VTKProperties& rhs)
{
  VTKProperties temp(rhs);
  Swap(temp);
  return *this;
}

void VTKProperties::Swap(VTKProperties& other)
{
  std::swap(m_lowerThreshold, other.m_lowerThreshold);
  std::swap(m_upperThreshold, other.m_upperThreshold);
  std::swap(m_mcThreshold, other.m_mcThreshold);
  std::swap(m_interpMode, other.m_interpMode);
  std::swap(m_kernelX, other.m_kernelX);
  std::swap(m_kernelY, other.m_kernelY);
  std::swap(m_kernelZ, other.m_kernelZ);
  std::swap(m_relaxationFactor, other.m_relaxationFactor);
  std::swap(m_iterations, other.m_iterations);
  std::swap(m_ambient, other.m_ambient);
  std::swap(m_diffuse, other.m_diffuse);
  std::swap(m_opacity, other.m_opacity);
  std::swap(m_specular, other.m_specular);
  std::swap(m_specularPower, other.m_specularPower);
  std::swap(m_featureAngle, other.m_featureAngle);
  std::swap(m_colorR, other.m_colorR);
  std::swap(m_colorG, other.m_colorG);
  std::swap(m_colorB, other.m_colorB);
  std::swap(m_stdDev[0], other.m_stdDev[0]);
  std::swap(m_stdDev[1], other.m_stdDev[1]);
  std::swap(m_stdDev[2], other.m_stdDev[2]);
  std::swap(m_radii[0], other.m_radii[0]);
  std::swap(m_radii[1], other.m_radii[1]);
  std::swap(m_radii[2], other.m_radii[2]);
}

template <typename T>
struct DeleteVTKObject:
  public std::unary_function<T*, void> {

  void operator()(T* ptr) const
  {
    ptr->Delete();
  }
};

struct VTKWidget::Implementation
{
  Implementation() :   
    m_gaussian(vtkImageGaussianSmooth::New()),
    m_smooth(vtkSmoothPolyDataFilter::New()),
    m_normals(vtkPolyDataNormals::New()),
    m_thresh(vtkImageThreshold::New()),
    m_close(vtkImageOpenClose3D::New()),
    m_brainActor(vtkActor::New()),
    m_sorter(vtkDepthSortPolyData::New()),
    m_surfMapper(vtkPolyDataMapper::New()),
    m_clipper(vtkClipPolyData::New()),
    m_mc(vtkMarchingCubes::New()),
    m_erode(vtkImageDilateErode3D::New())
  {};

  ~Implementation()
  {
    TRACKER("VTKWidget::Implementation::~Implementation");

    MESSAGE("destroying implementation elements");

    std::for_each(m_actors.begin(), m_actors.end(), 
		  DeleteVTKObject<vtkActor>());
    std::for_each(m_gaussians.begin(), m_gaussians.end(), 
		  DeleteVTKObject<vtkImageGaussianSmooth>());

    m_gaussian->Delete();
    m_smooth->Delete();
    m_normals->Delete();
    m_thresh->Delete();
    m_close->Delete();
    m_brainActor->Delete();
    m_sorter->Delete();
    m_surfMapper->Delete();
    m_clipper->Delete();
    m_mc->Delete();
    m_erode->Delete();
  }

  vtkImageGaussianSmooth *m_gaussian;
  vtkSmoothPolyDataFilter *m_smooth;
  vtkPolyDataNormals *m_normals;
  vtkImageThreshold *m_thresh;
  vtkImageOpenClose3D *m_close;
  vtkActor *m_brainActor;
  vtkDepthSortPolyData *m_sorter;
  vtkPolyDataMapper *m_surfMapper;
  vtkClipPolyData *m_clipper;
  vtkMarchingCubes *m_mc;
  vtkImageDilateErode3D *m_erode;

  std::vector<vtkActor *> m_actors;
  std::vector<MetaImage::Handle> m_metaImages;
  std::vector<vtkImageGaussianSmooth *> m_gaussians;

  VTKProperties m_props;
};

VTKWidget::VTKWidget(QWidget *parent, 
		     ImageGroup::Handle i,
		     OverlayList::Handle ol, 
		     Cursor::Handle c) :  
  ImageWidget(parent,i,ol,c), m_image(i), m_impl(new Implementation())
{
  TRACKER("VTKWidget::VTKWidget");

  QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );

  m_vtkwidget = new QVTKWidget( this, "vtkWidget" );
  setCentralWidget(m_vtkwidget);

  QToolBar *tb = new QToolBar(this);
  VTKToolbar *vt = new VTKToolbar(tb);

  addDockWindow(tb, tr("VTK Rendering Tools"), Top, FALSE);  

  connect(vt->m_printButton,   SIGNAL(clicked()), SLOT(print()));
  connect(vt->m_optionsButton, SIGNAL(clicked()), SLOT(options()));

  vtkRenderer *ren = vtkRenderer::New();
  ren->SetBackground(0.4, 0.4, 0.4);
  m_vtkwidget->GetRenderWindow()->AddRenderer(ren);

  vtkImageData *baseImage = vtkImageData::New();

  baseImage->SetDimensions(m_image->inqX(), m_image->inqY(), m_image->inqZ());
  ImageInfo::Handle miInfo = m_image->getMainImage()->getInfo();
  baseImage->SetSpacing(miInfo->inqXDim(), miInfo->inqYDim(), miInfo->inqZDim());
  baseImage->SetScalarTypeToFloat();
  baseImage->SetNumberOfScalarComponents(1);
  baseImage->AllocateScalars();

  float *ptr = (float *) baseImage->GetScalarPointer();
  for(unsigned int z = 0; z < m_image->inqZ(); ++z)
    for(unsigned int y = 0; y < m_image->inqY(); ++y)
      for(unsigned int x = 0; x < m_image->inqX(); ++x)
	{
	  *ptr++ = 
	    m_image->getMainImage()->getVolume(0)->value(m_image->inqX() - 1 - x, y, z);
	  // NB. X inversion to view 3D in neurological convention since data is in radiological order
	}
  
  m_impl->m_thresh->ThresholdBetween(m_impl->m_props.inqLowerThreshold(), m_impl->m_props.inqUpperThreshold());
  m_impl->m_thresh->SetOutputScalarTypeToFloat();
  m_impl->m_thresh->ReplaceOutOn();
  m_impl->m_thresh->ReplaceInOn();
  m_impl->m_thresh->SetOutValue(0);
  m_impl->m_thresh->SetInValue(1);
  m_impl->m_thresh->SetInput(baseImage);//(m_impl->m_gaussian->GetOutput());
  m_impl->m_thresh->UpdateInformation();

  int kx(0), ky(0), kz(0);
  m_impl->m_props.inqKernel(kx, ky, kz);
//   m_impl->m_close->SetInput(m_impl->m_thresh->GetOutput());
//   m_impl->m_close->SetKernelSize(kx, ky, kz);
//   m_impl->m_close->SetOpenValue(1.0);
//   m_impl->m_close->SetCloseValue(0.0);

//   m_impl->m_erode->SetInput(m_impl->m_thresh->GetOutput());
//   m_impl->m_erode->SetKernelSize(kx, ky, kz);
//   m_impl->m_erode->SetErodeValue(1.0);
//   m_impl->m_erode->SetDilateValue(0.0);
	  
  float sdX, sdY, sdZ;
  float rX, rY, rZ;
  m_impl->m_props.inqStdDev(sdX, sdY, sdZ);
  m_impl->m_props.inqRadii(rX, rY, rZ);
  m_impl->m_gaussian->SetInput(m_impl->m_thresh->GetOutput());
  m_impl->m_gaussian->SetDimensionality(3);
  m_impl->m_gaussian->SetStandardDeviations(sdX, sdY, sdZ);
  m_impl->m_gaussian->SetRadiusFactors(rX, rY, rZ);

  m_impl->m_mc->SetValue(0, m_impl->m_props.inqMcThreshold());
  m_impl->m_mc->SetInput(m_impl->m_gaussian->GetOutput());

  vtkPolyDataConnectivityFilter *connect = vtkPolyDataConnectivityFilter::New();
  connect->SetInput(m_impl->m_mc->GetOutput());
  connect->SetExtractionModeToLargestRegion();

  vtkDecimatePro *decim = vtkDecimatePro::New();
  decim->SetInput(connect->GetOutput());
  //decim->SetTargetReduction(10000);
  decim->PreserveTopologyOn();
  decim->SplittingOn();
  //decim->BoundaryVertexDeletionOn();
  //decim->SetMaximumError(VTK_LARGE_FLOAT);

  m_impl->m_normals->SetInput(decim->GetOutput());
  m_impl->m_normals->SplittingOn();
  m_impl->m_normals->ConsistencyOn();
  m_impl->m_normals->NonManifoldTraversalOff();
  m_impl->m_normals->SetFeatureAngle(m_impl->m_props.inqFeatureAngle());

  m_impl->m_smooth->SetInput(m_impl->m_normals->GetOutput());
  m_impl->m_smooth->SetNumberOfIterations(m_impl->m_props.inqIterations());
  m_impl->m_smooth->SetRelaxationFactor(m_impl->m_props.inqRelaxationFactor());

  VTKLookUpTable *lut = new VTKLookUpTable();

  vtkImageMapToColors *baseLayerRGBA = vtkImageMapToColors::New();
  baseLayerRGBA->SetOutputFormatToRGBA();
  baseLayerRGBA->SetInput(baseImage);
  baseLayerRGBA->SetLookupTable(VTKLookUpTable::GreyScale());
  float min(ol->getMainMetaImage()->getDs()->inqBriCon()->inqMin());
  float max(ol->getMainMetaImage()->getDs()->inqBriCon()->inqMax());
  baseLayerRGBA->GetLookupTable()->SetRange(min, max);

//   {
//     vtkDataSetWriter *dw = vtkDataSetWriter::New();
//     dw->SetInput(baseLayerRGBA->GetOutput());
//     dw->SetFileName("/tmp/baseLayerRGBA");
//     dw->Write();
//     dw->Delete();
//   }

  vtkImageBlend *blend = vtkImageBlend::New();
  blend->SetBlendModeToNormal();
  blend->AddInput(baseLayerRGBA->GetOutput());

  unsigned int count(1);
  MetaImageListIt it = ol->begin();
  ++it;
  for(; it != ol->end(); ++it)
    {
      // Process the overlays for possible rendering/blending along
      // with the main image surface.
      MetaImage::Handle mi = (*it);

      if(mi->inqVisibility()) {

	Image::Handle im(mi->getImage());
	ImageInfo::Handle info(mi->getInfo());
	BriCon::Handle bc(mi->getDs()->inqBriCon());

	vtkImageData *layer = vtkImageData::New();
	layer->SetDimensions(info->inqX(), info->inqY(), info->inqZ());
	layer->SetSpacing(info->inqXDim(), info->inqYDim(), info->inqZDim());
	layer->SetScalarTypeToFloat();
	layer->SetNumberOfScalarComponents(1);
	layer->AllocateScalars();
	
	vtkPoints *points = vtkPoints::New();
	vtkFloatArray *tensors = vtkFloatArray::New();
	tensors->SetNumberOfComponents(9);

	unsigned int offset(0);

	if(info->isMaskImage()) {
	  // Examine the image to determine what labels exist
	  // then loop over them to add each as new surface.
	  typedef std::map<float, bool> LabelMap;
	  LabelMap labels;

	  float *ptr = (float *)layer->GetScalarPointer();
	  for(unsigned int z = 0; z < info->inqZ(); ++z)
	    for(unsigned int y = 0; y < info->inqY(); ++y)
	      for(unsigned int x = 0; x < info->inqX(); ++x)
		{
		  float v = im->getVolume(0)->value(info->inqX() - 1 - x, y, z); 
		  *ptr++ = v;
		  if(v != 0) labels[v] = true;
		}
	  
	  // Convert image to RGBA and blend it into the viewed volume
	  MESSAGE("Blending in mask image");
	  vtkImageMapToColors *layerRGBA = vtkImageMapToColors::New();
	  layerRGBA->SetOutputFormatToRGBA();
	  layerRGBA->SetInput(layer);
	  LookUpTable::Handle lh(mi->getDs()->inqLookUpTable());
	  vtkLookupTable *thisLut;
	  if(lh)
	    thisLut = lut->GetLookupTableByName(lh->inqLutName());
	  else
	    thisLut = VTKLookUpTable::GreyScale();
	  thisLut->SetRange(bc->inqMin(), bc->inqMax());
	  MESSAGE(QString("Range: %1 to %2").arg(bc->inqMin()).arg(bc->inqMax()));
	  layerRGBA->SetLookupTable(thisLut);
	  blend->AddInput(layerRGBA->GetOutput());
	  blend->SetOpacity(count, 1.0); count++;

	  // Iterate over the labels creating a surface for each one.
	  for(LabelMap::iterator pos = labels.begin(); pos != labels.end(); ++pos)
	    {
	      float label = pos->first;

	      MESSAGE(QString("...rendering surface for mask label = %1").arg(label));

// 	      vtkImageThreshold *th = vtkImageThreshold::New();
// 	      th->SetInput(layer);
// 	      th->ThresholdBetween((label-0.5), (label+0.5));
// 	      th->ReplaceOutOn();
// 	      th->ReplaceInOn();
// 	      th->SetOutputScalarTypeToUnsignedChar();
// 	      th->SetOutValue(0.0);
// 	      th->SetInValue(label);
// 	      th->UpdateInformation();

	      vtkImageGaussianSmooth *gaussian = vtkImageGaussianSmooth::New();	      
	      gaussian->SetInput(layer);
	      gaussian->SetDimensionality(3);
	      gaussian->SetStandardDeviations(sdX, sdY, sdZ);
	      gaussian->SetRadiusFactors(rX, rY, rZ);
	      m_impl->m_gaussians.push_back(gaussian);

	      vtkImageThreshold *th = vtkImageThreshold::New();
	      th->SetInput(gaussian->GetOutput());
	      th->ThresholdBetween((label-0.5), (label+0.5));
	      th->ReplaceOutOn();
	      th->ReplaceInOn();
	      th->SetOutputScalarTypeToUnsignedChar();
	      th->SetOutValue(0.0);
	      th->SetInValue(label);
	      th->UpdateInformation();

	      vtkMarchingCubes *mclayer = vtkMarchingCubes::New();
	      mclayer->SetInput(gaussian->GetOutput());
	      mclayer->SetValue(0,label-.5);
	  
	      vtkPolyDataNormals *normals = vtkPolyDataNormals::New();
	      normals->SetInput(mclayer->GetOutput());
	      normals->SplittingOn();
	      normals->ConsistencyOn();
	      normals->NonManifoldTraversalOff();
	      normals->SetFeatureAngle(m_impl->m_props.inqFeatureAngle());

	      vtkPolyDataMapper *layerMapper = vtkPolyDataMapper::New();
	      layerMapper->SetInput(normals->GetOutput());
	      layerMapper->SetColorModeToMapScalars();
 	      layerMapper->UseLookupTableScalarRangeOn();
 	      layerMapper->SetLookupTable(thisLut);
	      layerMapper->ScalarVisibilityOn();

	      vtkActor *layerActor = vtkActor::New();
	      layerActor->SetMapper(layerMapper);
	      layerActor->GetProperty()->SetOpacity(0.4);

	      ren->AddViewProp(layerActor);
	  
	      m_impl->m_actors.push_back(layerActor);
	      m_impl->m_metaImages.push_back(mi);
	    }
	}
	if(info->isStatImage()) {
	  // Accumulate this layer into composite "blend"
	  // image for rendering onto cut surfaces
	  MESSAGE("Blending in stat image");
	  float *ptr = (float *)layer->GetScalarPointer();
	  for(unsigned int z = 0; z < info->inqZ(); ++z)
	    for(unsigned int y = 0; y < info->inqY(); ++y)
	      for(unsigned int x = 0; x < info->inqX(); ++x)
		*ptr++ = im->getVolume(0)->value(info->inqX() - 1 - x, y, z);
	  // Convert image to RGBA
	  vtkImageMapToColors *layerRGBA = vtkImageMapToColors::New();
	  layerRGBA->SetOutputFormatToRGBA();
	  layerRGBA->SetInput(layer);
	  LookUpTable::Handle lh(mi->getDs()->inqLookUpTable());
	  vtkLookupTable *layerLut, *surfLut;
	  if(lh)
	    {
	      layerLut = lut->GetLookupTableByName(lh->inqLutName());
	      surfLut  = lut->GetLookupTableByName(lh->inqLutName());
	    }
	  else
	    {
	      layerLut = VTKLookUpTable::GreyScale();
	      surfLut  = VTKLookUpTable::GreyScale();
	    }
	  layerLut->SetRange(bc->inqMin(), bc->inqMax());
	  layerRGBA->SetLookupTable(layerLut);
	  
	  blend->AddInput(layerRGBA->GetOutput());
	  blend->SetOpacity(count, 1.0); count++;
	  
	  // And IsoSurface
// 	  vtkImageThreshold *th = vtkImageThreshold::New();
// 	  th->SetInput(layer);
// 	  th->ThresholdByUpper(3);
// 	  th->ReplaceOutOn();
// 	  th->ReplaceInOn();
// 	  th->SetOutputScalarTypeToFloat();
// 	  th->SetOutValue(0.0);
// 	  th->SetInValue(1.0);
// 	  th->UpdateInformation();

// 	  vtkImageDilateErode3D *mcdil = vtkImageDilateErode3D::New();
// 	  mcdil->SetInput(th->GetOutput());
// 	  mcdil->SetKernelSize(3,3,3);
// 	  mcdil->SetDilateValue(1);
// 	  mcdil->SetErodeValue(0);
	  
	  vtkMarchingCubes *mclayer = vtkMarchingCubes::New();
	  mclayer->SetInput(layer);
	  mclayer->SetValue(0,bc->inqMin());
// 	  mclayer->SetValue(1,bc->inqMax());

	  vtkSmoothPolyDataFilter *smoothLayer = vtkSmoothPolyDataFilter::New();
	  smoothLayer->SetInput(mclayer->GetOutput());
	  smoothLayer->SetNumberOfIterations(5);
	  smoothLayer->SetFeatureEdgeSmoothing(true);
	  smoothLayer->SetBoundarySmoothing(true);
	  smoothLayer->SetRelaxationFactor(0.5);

//   {
//     vtkPolyDataWriter *dw = vtkPolyDataWriter::New();
//     dw->SetInput(mclayer->GetOutput());
//     dw->SetFileName("/tmp/mclayer");
//     dw->Write();
//     dw->Delete();
//   }
	  vtkPolyDataMapper *layerMapper = vtkPolyDataMapper::New();
	  layerMapper->SetInput(smoothLayer->GetOutput());
	  //layerMapper->SetScalarRange(1,2);
	  layerMapper->SetColorModeToMapScalars();
	  layerMapper->UseLookupTableScalarRangeOn();
	  surfLut->SetRange(0,bc->inqMin()*64);//bc->inqMin(), bc->inqMax());
	  layerMapper->SetLookupTable(surfLut);
	  layerMapper->ScalarVisibilityOn();

	  vtkActor *layerActor = vtkActor::New();
	  layerActor->SetMapper(layerMapper);
	  //bg 	  layerActor->GetProperty()->SetColor(surfLut->GetTableValue(16000));
	  layerActor->GetProperty()->SetOpacity(0.4);

	  ren->AddViewProp(layerActor);
	  
	  m_impl->m_actors.push_back(layerActor);
	  m_impl->m_metaImages.push_back(mi);
	  
	} else if(info->isDtiImage()) {
	  // Create glyphs for each tensor
	  // in the data set.
	  MESSAGE("Showing tensor image");

	  float *ptr = (float *)layer->GetScalarPointer();
	  for(unsigned int z = 0; z < info->inqZ(); ++z)
	    for(unsigned int y = 0; y < info->inqY(); ++y)
	      for(unsigned int x = 0; x < info->inqX(); ++x)
		{
		  float vx, vy, vz, mmx, mmy, mmz;

		  vx = im->getVolume(0)->value(x, y, z);
		  vy = im->getVolume(1)->value(x, y, z);
		  vz = im->getVolume(2)->value(x, y, z);

		  if((vx != 0) && (vy != 0) && (vz != 0)) {
		    FslGetMMCoord(info->inqStdMat(), 
				  x, y, z, &mmx, &mmy, &mmz);
		  
		    float tensor[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
		    for(unsigned short j = 0; j < 9; ++j)
		      tensor[j] = im->getVolume(j)->value(x, y, z);

		    points->InsertPoint(offset, mmx, mmy, mmz);
		    tensors->InsertTuple(offset, tensor);
		    ++offset;
		  }
		}

	  vtkSphereSource *sphere = vtkSphereSource::New();
	  sphere->SetCenter(0.0, 0.0, 0.0);
	  sphere->SetRadius(1.0);
	  sphere->SetThetaResolution(20);
	  sphere->SetPhiResolution(20);

	  vtkPolyData *tensorData = vtkPolyData::New();
	  tensorData->SetPoints(points); points->Delete();
	  tensorData->GetPointData()->SetTensors(tensors); tensors->Delete();

	  vtkTensorGlyph *glyph = vtkTensorGlyph::New();
	  glyph->ExtractEigenvaluesOff();
// 	  glyph->ThreeGlyphsOn();
	  glyph->SetInput(tensorData);
	  glyph->SetSource(sphere->GetOutput());

	  vtkPolyDataNormals *normals = vtkPolyDataNormals::New();
	  normals->SetInput(glyph->GetOutput());
	  vtkPolyDataMapper *ellipseMapper = vtkPolyDataMapper::New();
	  ellipseMapper->SetInput(normals->GetOutput());
	  
	  LookUpTable::Handle lh(mi->getDs()->inqLookUpTable());
	  vtkLookupTable *thisLut;
	  if(lh)
	    thisLut = lut->GetLookupTableByName(lh->inqLutName());
	  else
	    thisLut = VTKLookUpTable::GreyScale();
	       

	  vtkActor *ellipseActor = vtkActor::New();
	  ellipseActor->SetMapper(ellipseMapper);
	  ellipseActor->GetProperty()->SetOpacity(0.8);
	  ellipseActor->GetProperty()->SetColor(thisLut->GetTableValue(16000));

	  ren->AddViewProp(ellipseActor);
	  m_impl->m_actors.push_back(ellipseActor);
	  m_impl->m_metaImages.push_back(mi);
	}
      }
    }

//   vtkPiecewiseFunction *tfun = vtkPiecewiseFunction::New();
//   tfun->AddPoint(500,  0.0);
//   tfun->AddPoint(32000,0.8);

//   vtkColorTransferFunction *ctfun = vtkColorTransferFunction::New();
//   ctfun->AddRGBPoint(500.0,   0.0, 0.0, 0.0);
//   ctfun->AddRGBPoint(32000.0, 1.0, 1.0, 1.0);
  
//   vtkVolumeProperty *vprop = vtkVolumeProperty::New();
//   vprop->SetColor(ctfun);
//   vprop->SetScalarOpacity(tfun);

//   vtkVolumeRayCastCompositeFunction *cf = 
//     vtkVolumeRayCastCompositeFunction::New();

//   // Mapper
//   vtkVolumeRayCastMapper *mapper = vtkVolumeRayCastMapper::New();
//   mapper->SetVolumeRayCastFunction(cf);
//   mapper->SetInput(baseImage);

//   // Actor in scene
//   vtkVolume *vol = vtkVolume::New();
//   vol->SetMapper(mapper);
//   vol->SetProperty(vprop);  
//   ren->AddViewProp(vol);

  vtkBoxWidget *box = vtkBoxWidget::New();
  box->SetKeyPressActivationValue('b');
  box->SetInteractor(ren->GetRenderWindow()->GetInteractor());
  box->PlaceWidget(10,200,10,200,10,200);
  vtkPlanes *planes = vtkPlanes::New();
  box->GetPlanes(planes);
  vtkPolyData *pd = vtkPolyData::New();
  box->GetPolyData(pd);

  m_impl->m_clipper->SetInput(m_impl->m_smooth->GetOutput());
  m_impl->m_clipper->SetClipFunction(planes);

  m_impl->m_sorter->SetInput(m_impl->m_clipper->GetOutput());
  m_impl->m_sorter->SetCamera(ren->GetActiveCamera());
  m_impl->m_sorter->SortScalarsOn();
  m_impl->m_sorter->SetDirectionToBackToFront();

  if(m_impl->m_props.inqOpacity() != 1.0)
    m_impl->m_surfMapper->SetInput(m_impl->m_sorter->GetOutput());
  else
    m_impl->m_surfMapper->SetInput(m_impl->m_clipper->GetOutput());
  m_impl->m_surfMapper->ScalarVisibilityOff();

//   vtkTransform *boxXForm = vtkTransform::New();
//   box->GetTransform(boxXForm);
//   xform->SetTransform(boxXForm);

  vtkTransform *boxXForm = vtkTransform::New();
  box->GetTransform(boxXForm);

  ClippingBoxCallback *callback = ClippingBoxCallback::New();
  callback->SetPlanes(planes);
  callback->SetTransform(boxXForm);
  box->AddObserver(vtkCommand::EndInteractionEvent, callback);

  InteractorCallback *icallback = InteractorCallback::New();
  icallback->SetClipper(m_impl->m_clipper);
  icallback->SetMapper(m_impl->m_surfMapper);
  icallback->SetSorter(m_impl->m_sorter);
  icallback->SetProperties(&(m_impl->m_props));
  m_vtkwidget->GetInteractor()->AddObserver(vtkCommand::StartInteractionEvent, icallback);
  m_vtkwidget->GetInteractor()->AddObserver(vtkCommand::EndInteractionEvent, icallback);

  m_impl->m_brainActor->SetMapper(m_impl->m_surfMapper);
  m_impl->m_brainActor->GetProperty()->SetInterpolation(m_impl->m_props.inqInterpMode());
  m_impl->m_brainActor->GetProperty()->SetAmbient(m_impl->m_props.inqAmbient());
  m_impl->m_brainActor->GetProperty()->SetDiffuse(m_impl->m_props.inqDiffuse());
  m_impl->m_brainActor->GetProperty()->SetOpacity(m_impl->m_props.inqOpacity());
  m_impl->m_brainActor->GetProperty()->SetSpecular(m_impl->m_props.inqSpecular());
  m_impl->m_brainActor->GetProperty()->SetSpecularPower(m_impl->m_props.inqSpecularPower());
  float cr, cg, cb;
  m_impl->m_props.inqColor(cr, cg, cb);
  m_impl->m_brainActor->GetProperty()->SetColor(cr, cg, cb);


  // Add Actor to renderer
  ren->AddViewProp(m_impl->m_brainActor);

  unsigned short verts[][3] = { {1, 5, 0}, 
				{0, 4, 3}, {1, 2, 0}, {2, 6, 3},
				{5, 6, 4}, {1, 2, 5}
  };

//   {
//     vtkDataSetWriter *dw = vtkDataSetWriter::New();
//     dw->SetInput(blend->GetOutput());
//     dw->SetFileName("/tmp/blend2");
//     dw->Write();
//     dw->Delete();
//   }

//   vtkPolyDataToImageStencil *pdts = vtkPolyDataToImageStencil::New();
//   pdts->SetInput(m_impl->m_smooth->GetOutput());  
//   vtkImageStencil *stencil = vtkImageStencil::New();
//   stencil->SetInput(blend->GetOutput());
//   stencil->SetStencil(pdts->GetOutput());
//   stencil->ReverseStencilOff();
//   stencil->SetBackgroundValue(0);

  MESSAGE("Adding clipped planes");
  for(unsigned short p = 0; p < 6; ++p)
    {
      vtkPlaneSource *ps = vtkPlaneSource::New();

      ps->SetOrigin(pd->GetPoints()->GetPoint(verts[p][0]));
      ps->SetPoint1(pd->GetPoints()->GetPoint(verts[p][1]));
      ps->SetPoint2(pd->GetPoints()->GetPoint(verts[p][2]));

      ps->SetXResolution(100);
      ps->SetYResolution(100);
  
      vtkTransformPolyDataFilter *xform = vtkTransformPolyDataFilter::New();
      xform->SetInput(ps->GetOutput());
      xform->SetTransform(boxXForm);

      vtkImplicitDataSet *clipFun = vtkImplicitDataSet::New();
      clipFun->SetDataSet(m_impl->m_thresh->GetOutput());
//      clipFun->SetDataSet(stencil->GetOutput());
      clipFun->SetOutValue(0.0);

      vtkClipPolyData *clipPoly = vtkClipPolyData::New();
      clipPoly->SetClipFunction(clipFun);
      clipPoly->SetInput(xform->GetOutput());

      vtkProbeFilter *probe = vtkProbeFilter::New();
      probe->SetInput(clipPoly->GetOutput());
      probe->SetSource(blend->GetOutput());

//       vtkPolyDataWriter *pdw = vtkPolyDataWriter::New();
//       pdw->SetInput(probe->GetPolyDataOutput());
//       pdw->SetFileName("/tmp/Goop");
//       pdw->Write();

      vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
      mapper->SetInput(probe->GetPolyDataOutput());
//       mapper->SetInput(clipPoly->GetOutput());

      vtkActor *cutActor = vtkActor::New();
      cutActor->SetMapper(mapper);
      cutActor->GetProperty()->SetAmbient(0.7);
      cutActor->GetProperty()->SetDiffuse(0.3);
      cutActor->GetProperty()->SetInterpolationToGouraud();
      cutActor->GetProperty()->SetOpacity(1);

      ren->AddViewProp(cutActor);
    }

  vtkImagePlaneWidget *planex = vtkImagePlaneWidget::New();
  planex->SetInput(m_impl->m_gaussian->GetOutput());
  planex->SetPlaneOrientationToXAxes();
  planex->GetColorMap()->GetLookupTable()->SetAlpha(0.6);
  planex->SetKeyPressActivationValue('x');
  planex->GetTexturePlaneProperty()->SetOpacity(1);
  planex->DisplayTextOn();
  vtkImagePlaneWidget *planey = vtkImagePlaneWidget::New();
  planey->SetInput(m_impl->m_gaussian->GetOutput());
  planey->SetPlaneOrientationToYAxes();
  planey->GetColorMap()->GetLookupTable()->SetAlpha(0.6);
  planey->SetKeyPressActivationValue('y');
  planey->GetTexturePlaneProperty()->SetOpacity(1);
  planey->DisplayTextOn();
  vtkImagePlaneWidget *planez = vtkImagePlaneWidget::New();
  planez->SetInput(m_impl->m_gaussian->GetOutput());
  planez->SetPlaneOrientationToZAxes();
  planez->GetColorMap()->GetLookupTable()->SetAlpha(0.6);
  planez->SetKeyPressActivationValue('z');
  planez->GetTexturePlaneProperty()->SetOpacity(1);
  planez->DisplayTextOn();

  planey->SetLookupTable(planex->GetLookupTable());
  planez->SetLookupTable(planex->GetLookupTable());

  planex->SetInteractor(ren->GetRenderWindow()->GetInteractor());
  planey->SetInteractor(ren->GetRenderWindow()->GetInteractor());
  planez->SetInteractor(ren->GetRenderWindow()->GetInteractor());

  MESSAGE("Rendering");
  
  // Reset camera
  ren->ResetCamera();
  m_vtkwidget->GetRenderWindow()->Render();
  QApplication::restoreOverrideCursor();
}

VTKWidget::~VTKWidget()
{
  TRACKER("VTKWidget::~VTKWidget");
}

void VTKWidget::update(const Cursor::Handle c)
{
  TRACKER("VTKWidget::update(const Cursor::Handle)");
}

void VTKWidget::update(const OverlayList *ol, OverlayListMsg msg)
{
  TRACKER("VTKWidget::update(const OverlayList *ol, OverlayListMsg msg)");
  
  unsigned short count(0);

  for(std::vector<vtkActor *>::iterator it = m_impl->m_actors.begin(); 
      it != m_impl->m_actors.end(); ++it)
    {
      vtkActor *thisActor = m_impl->m_actors.at(count);
      MetaImage::Handle mi = m_impl->m_metaImages.at(count);

      thisActor->GetProperty()->SetOpacity(mi->inqTransparency());
      ++count; 
      QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );
      m_vtkwidget->GetRenderWindow()->Render();
      QApplication::restoreOverrideCursor();
    }
}

void VTKWidget::print() 
{
  QString fn = QFileDialog::getSaveFileName("screenshot.tiff", 
					    "TIFF files (*.tif, *.tiff)", this,
					    "Screenshot dialog",
					    "Select a filename for saving");

  if(!fn.isNull()) {
    vtkWindowToImageFilter *w2i = vtkWindowToImageFilter::New();
    vtkTIFFWriter *writer = vtkTIFFWriter::New();
    
    w2i->SetInput(m_vtkwidget->GetRenderWindow());
    w2i->Update();
    writer->SetInput(w2i->GetOutput());
    writer->SetFileName((const char *)fn);
    QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );
    m_vtkwidget->GetRenderWindow()->Render();
    QApplication::restoreOverrideCursor();
    writer->Write();

    writer->Delete();
    w2i->Delete();
  }
}

struct SetGaussianParams
{
  SetGaussianParams(float sdx, float sdy, float sdz,
		    float  rx, float  ry, float  rz):
    m_sdx(sdx), m_sdy(sdy), m_sdz(sdz),
     m_rx( rx),  m_ry( ry),  m_rz( rz)
  {}

  void operator()(vtkImageGaussianSmooth *g) const
  {
    g->SetStandardDeviations(m_sdx, m_sdy, m_sdz);
    g->SetRadiusFactors(m_rx, m_ry, m_rz);
  }

  float m_sdx, m_sdy, m_sdz, m_rx, m_ry, m_rz;
};

void VTKWidget::options()
{
  VTKPropertyDialog optionsDialog(this, m_impl->m_props);
  
  if(optionsDialog.exec() == QDialog::Accepted)
    {
      m_impl->m_props = optionsDialog.getProperties();
      m_impl->m_thresh->ThresholdBetween(m_impl->m_props.inqLowerThreshold(), m_impl->m_props.inqUpperThreshold());
      float sdX, sdY, sdZ;
      m_impl->m_mc->SetValue(0, m_impl->m_props.inqMcThreshold());
      m_impl->m_props.inqStdDev(sdX, sdY, sdZ);
      m_impl->m_gaussian->SetStandardDeviations(sdX, sdY, sdZ);
      float rX, rY, rZ;
      m_impl->m_props.inqRadii(rX, rY, rZ);
      m_impl->m_gaussian->SetRadiusFactors(rX, rY, rZ);

      std::for_each(m_impl->m_gaussians.begin(), m_impl->m_gaussians.end(),
		    SetGaussianParams(sdX, sdY, sdZ, rX, rY, rZ));

      m_impl->m_brainActor->GetProperty()->SetInterpolation(m_impl->m_props.inqInterpMode());
      m_impl->m_brainActor->GetProperty()->SetAmbient(m_impl->m_props.inqAmbient());
      m_impl->m_brainActor->GetProperty()->SetDiffuse(m_impl->m_props.inqDiffuse());
      m_impl->m_brainActor->GetProperty()->SetOpacity(m_impl->m_props.inqOpacity());
      m_impl->m_brainActor->GetProperty()->SetSpecular(m_impl->m_props.inqSpecular());
      m_impl->m_brainActor->GetProperty()->SetSpecularPower(m_impl->m_props.inqSpecularPower());
      float cr, cg, cb;
      m_impl->m_props.inqColor(cr, cg, cb);
      m_impl->m_brainActor->GetProperty()->SetColor(cr, cg, cb);
      int x, y, z;
      m_impl->m_props.inqKernel(x, y, z);
      //m_impl->m_close->SetKernelSize(x, y, z);
      //m_impl->m_erode->SetKernelSize(x, y, z);
      m_impl->m_smooth->SetNumberOfIterations(m_impl->m_props.inqIterations());
      m_impl->m_smooth->SetRelaxationFactor(m_impl->m_props.inqRelaxationFactor());
      m_impl->m_normals->SetFeatureAngle(m_impl->m_props.inqFeatureAngle());
      if(m_impl->m_props.inqOpacity() != 1.0)
	m_impl->m_surfMapper->SetInput(m_impl->m_sorter->GetOutput());
      else
	m_impl->m_surfMapper->SetInput(m_impl->m_clipper->GetOutput());

      QApplication::setOverrideCursor( QCursor(Qt::WaitCursor) );
      m_vtkwidget->GetRenderWindow()->Render();
      QApplication::restoreOverrideCursor();
    }
}

