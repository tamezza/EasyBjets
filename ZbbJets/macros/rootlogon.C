// rootlogon.C - run automatically by ROOT when it starts in ZbbJets/macros.
// Adds ../../common to the include path, so the headers in common/macros find
// VariantConfig.h, CommonUtils.h, ... This is why the macros must be run from
// this folder.
{
    gInterpreter->AddIncludePath("../../common");
    gInterpreter->AddIncludePath(".");
}
