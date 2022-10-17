/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

// This file declares functions for setting handlers (initialization, error, warning).

#ifndef __HANDLER_H__
#define __HANDLER_H__

#include <cstdarg>
#include <iostream>

typedef	void (*ZErrorHandler)(const char* module, const char* fmt, va_list);
void ZWarning(const char* module, const char* fmt, ...);
void ZError(const char* module, const char* fmt, ...);
ZErrorHandler ZSetWarningHandler(ZErrorHandler handler);
ZErrorHandler ZSetErrorHandler(ZErrorHandler handler);

#ifdef _MSC_VER
void SetWindowsHandler();
#endif

#ifdef _WINDOWS
	#include <afxwin.h>         // MFC core and standard components
	#include <afxext.h>
	#include <afxcmn.h>
	#include "..\mfclib\ProgressBar.h"
#endif

class ZProgressIndicator
{
#ifdef _WINDOWS
	CProgressBar bar;
#endif		

	const char* m_txt;
	int			m_max;
	float		inc, step;
	bool		m_enable;
public:
	static bool	m_progress;

	ZProgressIndicator(const char* txt, int size, int max, bool enable=true) : m_txt(txt), m_max(max), m_enable(enable)
	{
		if(!m_progress) return;

		inc = float(size) / max;
		step = 0;

#ifdef _WINDOWS
		bar.Create(txt, 80, max);
#else
		if(m_enable) std::cerr << "    " << txt << '(' << max << ')';
#endif
	}
	
	void Step(int its)
	{
		if(!m_progress) return;
#ifdef _WINDOWS
		CString str;
		str.Format("%s: %d%% complete", m_txt, its*100/m_max);
		bar.SetText(str);
		bar.StepIt();
#else
		if(m_enable) 
		{
			if((step+=inc)>1) step--, std::cerr << '.';
		}
#endif		
	}
};

inline void EnableProgressIndicator(bool prog = true) { ZProgressIndicator::m_progress = prog; }

#endif
