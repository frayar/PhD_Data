#ifndef EXPORT_H
#define EXPORT_H

	/* PROTOTYPES */
	void ExportGraph(std::map<__int64, double>* pList, const char *filePath, int row, EXPORT method);
	void WriteAdjacencies(std::map<__int64, double>* pList, const char *filePath, int row);
	void WriteJSONResults(std::map<__int64, double>* pList, const char *filePath, int row);
	void WriteGEXFResults(std::map<__int64, double>* pList, const char *filePath, int row);
	void WriteTLPResults(std::map<__int64, double>* pList, const char *filePath, int row);

	void ExportGraph(list<DoubleResultPackage>* pList, const char *filePath, int row, EXPORT method);

#endif