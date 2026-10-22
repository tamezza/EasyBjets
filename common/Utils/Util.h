#ifndef UTIL_H
#define UTIL_H

/*
 *  A bunch of useful methods for ROOT data analysis
 *  Dag Gillberg
 */

typedef TString Str;
typedef std::vector<TString> StrV;
typedef std::vector<double> NumV;
typedef std::vector<int> VecI;
typedef std::vector<int> IntV;
typedef std::map<TString, TH1D*> HistMap;
typedef std::vector<TH1D*> HistV;
typedef std::vector<HistV> HistVV;
typedef std::map<Str, NumV> NumVMap;
typedef TGraphAsymmErrors TGA;

TCanvas *_can;
map<Str,TFile*> _fMap;
double _lm = 0.12, _rm = 0.02, _top = 0.05, _bot = 0.12;

void fatal(Str msg) {
  printf("FATAL:\n\n  %s\n\n", msg.Data());
  exit(EXIT_FAILURE);
}

bool fileExist(Str fn) {
  return !(gSystem->AccessPathName(fn.Data()));
}

TFile *openFile(Str fn) {
  TFile *f = TFile::Open(fn);
  if (f == NULL || f->IsZombie())
    fatal("Cannot open " + fn);
  return f;
}

TFile *getFile(Str fn) {
  if (_fMap[fn]==nullptr) _fMap[fn]=openFile(fn);
  return _fMap[fn];
}

TTree *getTree(TFile *f, Str tn) {
  TTree *t = (TTree*)f->Get(tn);
  if (t == nullptr) fatal("Cannot access " + tn + " in " + f->GetName());
  return t;
}

TH1D *getHisto(TFile *f, Str hn) {
  TH1D *h = (TH1D*) f->Get(hn);
  if (h == nullptr) fatal("Cannot access " + hn + " in " + f->GetName());
  return h;
}

TH2D *getHisto2D(TFile *f, Str hn) {
  TH2D *h = (TH2D*) f->Get(hn);
  if (h == nullptr) fatal("Cannot access " + hn + " in " + f->GetName());
  return h;
}

TH2D *getHisto2D(Str fn, Str hn) {
  return getHisto2D(getFile(fn),hn);
}

StrV merge(StrV a, StrV b) {
  StrV v(a); for (auto e:b) v.push_back(e); return v;
}

TEnv *openSettingsFile(Str fileName) {
  if (fileName == "")
    fatal("No config file name specified. Cannot open file!");
  TEnv *settings = new TEnv();
  int status = settings->ReadFile(fileName.Data(), EEnvLevel(0));
  if (status != 0)
    fatal(Form("Cannot read file %s", fileName.Data()));
  return settings;
}

TH1D *getHisto(Str fn, Str hn) {
  return getHisto(getFile(fn), hn);
}

TTree *getTree(Str fn, Str tn) {
  return getTree(getFile(fn), tn);
}

StrV vectorize(Str str, Str sep=" ") {
  StrV result;
  TObjArray *strings = str.Tokenize(sep.Data());
  if (strings->GetEntries() == 0)
    return result;
  TIter istr(strings);
  while (TObjString* os = (TObjString*) istr())
    if (os->GetString()[0] != '#')
      result.push_back(os->GetString());
    else
      break;
  return result;
}

IntV interize(Str str, Str sep=" ") {
  IntV result;
  TObjArray *strings = str.Tokenize(sep.Data());
  if (strings->GetEntries() == 0)
    return result;
  TIter istr(strings);
  while (TObjString* os = (TObjString*) istr())
    if (os->GetString()[0] != '#')
      result.push_back((os->GetString()).Atoi());
    else
      break;
  return result;
}

NumV numberize(Str str, Str sep=" ") {
  NumV result; StrV words=vectorize(str,sep);
  for (auto s:words) result.push_back(s.Atof());
  return result;
}

Str getStr(TEnv *settings, Str key) {
  Str val = settings->GetValue(key,"");
  if (val=="") fatal("No value found for TEnv key: "+key);
  return val;
}

Str getStr(TEnv *settings, Str key, Str def) {
  return settings->GetValue(key,def);
}

StrV getStrV(TEnv *settings, Str key, Str sep=" ") {
  return vectorize(getStr(settings,key),sep);
}

double getNum(TEnv *settings, Str key) {
  double num = settings->GetValue(key,-999.99);
  if (num==-999.99) fatal("No value found for TEnv key: "+key);
  return num;
}

double getNum(TEnv *settings, Str key, double def) {
  return settings->GetValue(key,def);
}

double getInt(TEnv *settings, Str key) {
  int num = settings->GetValue(key,-999);
  if (num==-999) fatal("No value found for TEnv key: "+key);
  return num;
}

// Set style similar to ATLAS style
void setStyle() {
  double tsize = 0.04;
  gStyle->SetTextSize(tsize); gStyle->SetTextFont(42);
  for (auto ax:{"x","y"}) {
    gStyle->SetLabelSize(tsize,ax); gStyle->SetTitleSize(tsize,ax);
    gStyle->SetLabelFont(42,ax); gStyle->SetTitleFont(42,ax);
    gStyle->SetTickLength(0.02,ax);
  }
  gStyle->SetPadLeftMargin(_lm); gStyle->SetPadRightMargin(_rm);
  gStyle->SetPadBottomMargin(_bot); gStyle->SetPadTopMargin(_top);

  TH1::SetDefaultSumw2(); gStyle->SetHistMinimumZero();
  gStyle->SetOptTitle(0); gStyle->SetOptStat(0); gStyle->SetOptFit(0);

  gStyle->SetPadTickX(1); gStyle->SetPadTickY(1); // tick marks on top and RHS

  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
}

void drawText(double x, double y, Str txt, int col=kBlack, double tsize=0.04, int align=11) {
  static TLatex *tex = new TLatex(); tex->SetTextAlign(align);
  tex->SetTextColor(col); tex->SetNDC(); tex->SetTextSize(tsize); tex->SetTextFont(42);
  tex->DrawLatex(x, y, txt);
}

void drawCText(double x, double y, Str txt, int col=kBlack, double tsize=0.04) { drawText(x,y,txt,col,tsize,23); }
void drawRText(double x, double y, Str txt, int col=kBlack, double tsize=0.04) { drawText(x,y,txt,col,tsize,33); }
void drawVertText(double x, double y, Str txt, int col=kBlack, double tsize=0.04) {
  static TLatex *tex = new TLatex(); tex->SetTextAlign(33); tex->SetTextColor(col); tex->SetNDC();
  tex->SetTextSize(tsize); tex->SetTextFont(42); tex->SetTextAngle(90); tex->DrawLatex(x, y, txt);
}

void drawLine(double x1, double y1, double x2, double y2, int col = kBlack, int ls = 1, int lw = 1) {
  static TLine *line = new TLine();
  line->SetLineColor(col); line->SetLineStyle(ls); line->SetLineWidth(lw);
  line->DrawLine(x1, y1, x2, y2);
}

void drawLineNDC(double x1, double y1, double x2, double y2, int col = kBlack, int ls = 1, int lw = 1) {
  static TLine *line = new TLine();
  line->SetLineColor(col); line->SetLineStyle(ls); line->SetLineWidth(lw);
  line->DrawLineNDC(x1,y1,x2,y2);
}

void formatHisto(TH1 *h, int col=kAzure, int lw=2, int ls=1, int ms=1, float msize=1) {
  h->SetLineColor(col); h->SetMarkerColor(col);
  h->SetLineWidth(lw); h->SetLineStyle(ls); h->SetMarkerStyle(ms); h->SetMarkerSize(msize);
  h->SetStats(0);
}

void formatHistos(vector<TH1D*> hs, int col=kAzure, int lw=2, int ls=1, int ms=1) {
  for (size_t i = 0; i < hs.size(); ++i)
    formatHisto(hs[i], col, lw, ls, ms);
}

int getCol(int style) {
  static VecI cols({kBlack,kAzure,kRed,kGreen+1,kViolet,kOrange+2});
  return cols[style%6];
}

int getMS(int style) {
  static VecI ms({20,22,23,21,25,29,27,26});
  return ms[style%8];
}

void drawHisto(TH1 *h, Str opt="", int col=kBlack, int lw=2, int ls=1, int ms=0) {
  formatHisto(h,col,lw,ls,ms); h->Draw(opt);
}

void drawHistoStyle(TH1 *h, Str opt, int style) {
  int col=getCol(style), lw=2, ms=getMS(style); double msize=1.0; // style=0
  formatHisto(h,col,lw,1,ms); h->SetMarkerSize(msize); h->Draw(opt); h->Draw(opt);
}

void drawHistos(HistV hv, Str opt="") {
  if (hv.size()==0) return;
  drawHisto(hv[0],opt,getCol(1));
  for (size_t i=1;i<hv.size();++i)
    drawHisto(hv[i],"same"+opt,getCol(1+i),2,i/6+1);
}

void drawLegend(double x, double y, StrV descs, double tsize=0.04) {
  size_t N = descs.size();
  for (size_t i=0;i<N;++i)
    drawText(x,y-tsize*1.25*i,descs[i],getCol(i+1),tsize);
}

TH1D *draw1D(TTree *t, Str var, Str cut, int Nbins, double min, double max,
             int col=kBlack, int lw=2, int ls=1, int ms=0, Str tit="") {
  static int hi=0; Str hname(Form("hist%d",++hi));
  TH1D *h = new TH1D(hname,";"+var,Nbins,min,max);
  t->Draw(var+">>"+hname,cut,"e");
  formatHisto(h,col,lw,ls,ms);
  if (tit!="") h->SetTitle(tit);
  return h;
}

TH1D *draw1D(TTree *t, Str var, Str cut, NumV bins,
             int col=kBlack, int lw=2, int ls=1, int ms=0, Str tit="") {
  static int hi=0; Str hname(Form("histVB%d",++hi));
  TH1D *h = new TH1D(hname,";"+var,-1+bins.size(),&bins[0]);
  t->Draw(var+">>"+hname,cut,"e");
  formatHisto(h,col,lw,ls,ms);
  if (tit!="") h->SetTitle(tit);
  return h;
}

ROOT::RDF::RResultPtr<::TH1D> getRDF1D(ROOT::RDataFrame *df, Str var, Str cut, 
                                       int Nbins, double min, double max, Str wght="weight") {
  static int hi=0; Str hname(Form("df%d",++hi));
  ROOT::RDF::TH1DModel model={hname,";"+var,Nbins,min,max};
  auto h = df->Histo1D(model,(string_view)var);
  if (cut!="") h = df->Filter((string_view)cut).Histo1D(model,(string_view)var);
  return h;
}


TH2D *draw2D(TTree *t, Str xvar, Str yvar, Str cut, int NbinsX, double xmin, double xmax,
             int NbinsY, double ymin, double ymax) {
  static int hi=0; Str hname(Form("hist2D_%d",++hi));
  auto h = new TH2D(hname,";"+xvar+";"+yvar,NbinsX,xmin,xmax,NbinsY,ymin,ymax);
  t->Draw(yvar+":"+xvar+">>"+hname,cut,"colz");
  return h;
}

void prepTopCanvas(bool logy=false) {
  _can->cd(); _can->SetLogy(logy);
  _can->SetBottomMargin(0.36); _can->SetTopMargin(_top);
  _can->SetLeftMargin(_lm); _can->SetRightMargin(_rm);
}

void prepSubPad(bool logy=false) {
  TPad *sub = new TPad("","",0,0,1,1.0,0,0,0); sub->SetLogy(logy);
  sub->SetFillStyle(0);
  sub->SetTopMargin(0.64); sub->SetBottomMargin(_bot);
  sub->SetLeftMargin(_lm); sub->SetRightMargin(_rm); sub->Draw(); sub->cd();
}

double getLumi(Str fn) {
  if      (fn.Contains("2016")) return 35660.06;
  else if (fn.Contains("2017")) return 44307.4;
  else if (fn.Contains("2018")) return 58450.1;
  else return 35660.06+ 44307.4+58450.1;
  fatal("Can't figure out year for "+fn);
  return 0.0;
}

void printTime(bool printCur=true) {
  static bool first=true;
  static time_t start;
  if (first) { first=false; ::time(&start); }
  time_t aclock; ::time( &aclock );
  if (printCur) cout <<  "Current time: " << asctime( localtime( &aclock ) );
  cout  << " ( "<< ::difftime( aclock, start) <<" s elapsed )" << std::endl;
}

#endif
