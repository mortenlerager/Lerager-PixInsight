# Lerager PixInsight Scripts

AstroClean 2.8.7 and AlignRGB Preview 0.3, by Morten Lerager.

## PixInsight installation (test release, UNSIGNED)

In PixInsight 1.9.6: Resources > Updates > Manage Repositories > Add:

`https://raw.githubusercontent.com/mortenlerager/Lerager-PixInsight/main/updates/`

Then Resources > Updates > Check for Updates, apply and restart.

**Important:** This repository is currently unsigned. PixInsight may display a security/trust warning. Do not claim a signed installation until a PixInsight-registered signing identity is configured and tested.

The scripts are installed under `src/scripts/Lerager/` and appear in Script > Lerager.

## Releasing updates

Edit the script in `src/scripts/Lerager/`, update version strings as needed, run `python tools/build.py`, and upload the updated script and documentation, `updates/updates.xri`, and the newly built ZIP to GitHub. Each release must use a new package name (date alone is insufficient for multiple releases on the same day). Sign the XRI with the PixInsight-supported workflow once a registered identity is available.
