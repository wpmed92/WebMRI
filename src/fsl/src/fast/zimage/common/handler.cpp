/*  Part of FAST - FMRIB's Automated Segmentation Tool

    Yongyue Zhang, FMRIB Image Analysis Group

    Copyright (C) 1999-2003 University of Oxford  */

/*  CCOPYRIGHT */

#include <iostream>
#include <new>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "handler.h"
using namespace std;

#ifdef _MSC_VER
	static int new_failed(unsigned int)
	{
	#ifdef _WINDOWS
		AfxMessageBox("new: Out of Memory!");
	#else
		cerr << "new: Out of Memory!\n";
	#endif
		exit(0);
		return 0;
	}
#else
	static void new_failed()
	{
		cerr << "new: Out of Memory!\n";
		exit(0);
	}
#endif

static int init()
{
#ifdef _MSC_VER
	_set_new_handler(new_failed);
#else
	set_new_handler(new_failed);
#endif

	EnableProgressIndicator(true);

	return 0;
}

static int dummy=init();

bool ZProgressIndicator::m_progress = true;

#ifdef _MSC_VER

#include <windows.h>

static void WinWarningHandler(const char* module, const char* fmt, va_list ap)
{
	LPTSTR szTitle;
	LPTSTR szTmp;
	LPCTSTR szTitleText = "%s Warning";
	LPCTSTR szDefaultModule = "TIFFLIB";
	szTmp = (module == NULL) ? (LPTSTR)szDefaultModule : (LPTSTR)module;
	if ((szTitle = (LPTSTR)LocalAlloc(LMEM_FIXED, (lstrlen(szTmp) +
			lstrlen(szTitleText) + lstrlen(fmt) + 128)*sizeof(TCHAR))) == NULL)
		return;
	wsprintf(szTitle, szTitleText, szTmp);
	szTmp = szTitle + (lstrlen(szTitle)+2)*sizeof(TCHAR);
	vsprintf(szTmp, fmt, ap);
	MessageBox(GetFocus(), szTmp, szTitle, MB_OK | MB_ICONINFORMATION);
	LocalFree(szTitle);
}

static void WinErrorHandler(const char* module, const char* fmt, va_list ap)
{
	LPTSTR szTitle;
	LPTSTR szTmp;
	LPCTSTR szTitleText = "%s Error";
	LPCTSTR szDefaultModule = "TIFFLIB";
	szTmp = (module == NULL) ? (LPTSTR)szDefaultModule : (LPTSTR)module;
	if ((szTitle = (LPTSTR)LocalAlloc(LMEM_FIXED, (lstrlen(szTmp) +
			lstrlen(szTitleText) + lstrlen(fmt) + 128)*sizeof(TCHAR))) == NULL)
		return;
	wsprintf(szTitle, szTitleText, szTmp);
	szTmp = szTitle + (lstrlen(szTitle)+2)*sizeof(TCHAR);
	vsprintf(szTmp, fmt, ap);
	MessageBox(GetFocus(), szTmp, szTitle, MB_OK | MB_ICONEXCLAMATION);
	LocalFree(szTitle);
}

#endif

static void ConsoleWarningHandler(const char* module, const char* fmt, va_list ap)
{
	if (module != NULL)
		fprintf(stderr, "%s: ", module);
	fprintf(stderr, "Warning, ");
	vfprintf(stderr, fmt, ap);
	fprintf(stderr, "\n");
}

static void ConsoleErrorHandler(const char* module, const char* fmt, va_list ap)
{
	if (module != NULL)
		fprintf(stderr, "%s: ", module);
	vfprintf(stderr, fmt, ap);
	fprintf(stderr, "\n");
	exit(-1);
}

ZErrorHandler _ZWarningHandler = ConsoleWarningHandler;
ZErrorHandler _ZErrorHandler = ConsoleErrorHandler;

#ifdef _MSC_VER

void SetWindowsHandler()
{
	_ZWarningHandler = WinWarningHandler;
	_ZErrorHandler = WinErrorHandler;
}

#endif

ZErrorHandler ZSetWarningHandler(ZErrorHandler handler)
{
	ZErrorHandler prev = _ZWarningHandler;
	_ZWarningHandler = handler;
	return (prev);
}

void ZWarning(const char* module, const char* fmt, ...)
{
	if (_ZWarningHandler) 
	{
		va_list ap;
		va_start(ap, fmt);
		(*_ZWarningHandler)(module, fmt, ap);
		va_end(ap);
	}
}

ZErrorHandler ZSetErrorHandler(ZErrorHandler handler)
{
	ZErrorHandler prev = _ZErrorHandler;
	_ZErrorHandler = handler;
	return (prev);
}

void ZError(const char* module, const char* fmt, ...)
{
	if (_ZErrorHandler) 
	{
		va_list ap;
		va_start(ap, fmt);
		(*_ZErrorHandler)(module, fmt, ap);
		va_end(ap);
	}
}
