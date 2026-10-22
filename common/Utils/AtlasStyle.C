#include <TStyle.h>
#include "AtlasStyle.h"

void SetAtlasStyle() {
    static TStyle* atlasStyle = 0;
    if (!atlasStyle) {
        atlasStyle = AtlasStyle();
    }
    gROOT->SetStyle("ATLAS");
    gROOT->ForceStyle();
}

TStyle* AtlasStyle() {
    TStyle* atlasStyle = new TStyle("ATLAS", "Atlas style");

    // Set the paper & margin sizes
    atlasStyle->SetPaperSize(20, 26);

    // Set the margin sizes
    atlasStyle->SetPadTopMargin(0.05);
    atlasStyle->SetPadRightMargin(0.05);
    atlasStyle->SetPadBottomMargin(0.16);
    atlasStyle->SetPadLeftMargin(0.16);

    // Set title offsets (for axis label)
    atlasStyle->SetTitleXOffset(1.4);
    atlasStyle->SetTitleYOffset(1.4);

    // Use large fonts
    Int_t font = 42; // Helvetica
    Double_t tsize = 0.05;
    atlasStyle->SetTextFont(font);

    atlasStyle->SetTextSize(tsize);
    atlasStyle->SetLabelFont(font, "x");
    atlasStyle->SetTitleFont(font, "x");
    atlasStyle->SetLabelFont(font, "y");
    atlasStyle->SetTitleFont(font, "y");
    atlasStyle->SetLabelFont(font, "z");
    atlasStyle->SetTitleFont(font, "z");

    atlasStyle->SetLabelSize(tsize, "x");
    atlasStyle->SetTitleSize(tsize, "x");
    atlasStyle->SetLabelSize(tsize, "y");
    atlasStyle->SetTitleSize(tsize, "y");
    atlasStyle->SetLabelSize(tsize, "z");
    atlasStyle->SetTitleSize(tsize, "z");

    // Use bold lines and markers
    atlasStyle->SetMarkerStyle(20);
    atlasStyle->SetMarkerSize(1.2);
    atlasStyle->SetHistLineWidth(2.);
    atlasStyle->SetLineStyleString(2, "[12 12]"); // postscript dashes

    // Get rid of X error bars 
    // atlasStyle->SetErrorX(0.001);
    // Get rid of error bar caps
    atlasStyle->SetEndErrorSize(0.);

    // Do not display any of the standard histogram decorations
    atlasStyle->SetOptTitle(0);
    atlasStyle->SetOptStat(0);
    atlasStyle->SetOptFit(0);

    // Put tick marks on top and RHS of plots
    atlasStyle->SetPadTickX(1);
    atlasStyle->SetPadTickY(1);

    return atlasStyle;
}
