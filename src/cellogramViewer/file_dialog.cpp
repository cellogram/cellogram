///////////////////////////////////////////////////////////////////////////////
#include "file_dialog.h"
// -----------------------------------------------------------------------------
#include "tinyfiledialogs.h"
// -----------------------------------------------------------------------------
#include <memory>
#include <vector>
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

namespace FileDialog {

// -----------------------------------------------------------------------------

std::string openFileName(const std::string &defaultPath,
	const std::vector<std::string> &filters, const std::string &desc)
{
	int n = static_cast<int>(filters.size());
	std::vector<char const *> filterPatterns(n);
	for (int i = 0; i < n; ++i) {
		filterPatterns[i] = filters[i].c_str();
	}
	char const * select = tinyfd_openFileDialog("Open File",
		defaultPath.c_str(), n, filterPatterns.data(), desc.c_str(), 0);
	if (select == nullptr) {
		return "";
	} else {
		return std::string(select);
	}
}

std::string saveFileName(const std::string &defaultPath,
	const std::vector<std::string> &filters, const std::string &desc)
{
	int n = static_cast<int>(filters.size());
	std::vector<char const *> filterPatterns(n);
	for (int i = 0; i < n; ++i) {
		filterPatterns[i] = filters[i].c_str();
	}
	char const * select = tinyfd_saveFileDialog("Save File",
		defaultPath.c_str(), n, filterPatterns.data(), desc.c_str());
	if (select == nullptr) {
		return "";
	} else {
		return std::string(select);
	}
}

// -----------------------------------------------------------------------------

} // namespace FileDialog
