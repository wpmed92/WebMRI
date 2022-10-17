// testimage.cpp: implementation of the TestImage class.
//
//////////////////////////////////////////////////////////////////////

#include <assert.h>

#include "testimage.h"
#include "sum.h"

#include "../image.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

class TestImageImpl
{
public:
	TestImageImpl() {}
	~TestImageImpl() {}

	void test();
};

TestImage::TestImage()
{
	pImpl_ = new TestImageImpl;
}

TestImage::~TestImage()
{
	delete pImpl_;
}

void TestImage::run()
{
	pImpl_->test();
}

void TestImageImpl::test()
{
	Image::Handle i = Image::load("functional");

	assert(i->inqX() == 64);
	assert(i->inqY() == 64);
	assert(i->inqZ() == 21);
	assert(i->inqNumVolumes() == 180);

	Volume::Handle v = i->getVolume(0);

	Sum sum;
	v->accept(sum);
	
//	VolumeUS::Handle = v.downcast();

	float total = sum.total();
	int count = sum.count();

	assert(count == v->inqX() * v->inqY() * v->inqZ());

}