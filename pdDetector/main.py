import datetime
import logging
import glob
import os

import matlab.engine
import numpy as np
import fire

import iiic_detector.utils
import iiic_detector


def main(
    target_dir: str,
    output_dir: str,
    save_data: bool = False,
    save_png: bool = False,
    log: str = "WARNING",
):
    numeric_log_level = getattr(logging, log.upper())
    logging.basicConfig(level=numeric_log_level)

    date_str = datetime.date.today().strftime("%Y-%m-%d")

    analysis_dir = os.path.join(output_dir, f"analysis_{date_str}")
    image_dir = os.path.join(output_dir, "images")
    data_dir = os.path.join(analysis_dir, "data")

    iiic_detector.utils.make_dir_safe(analysis_dir)
    iiic_detector.utils.make_dir_safe(image_dir)
    iiic_detector.utils.make_dir_safe(data_dir)

    fnames = glob.glob(os.path.join(target_dir, "*.mat"))

    eng = matlab.engine.start_matlab()
    cutoff_pct = 1.0 / 3.0
    result_table = []

    for idx, fname in enumerate(fnames):
        logging.info(fname)

        # For testing
        if (
            fname
            != "../testdata/IIICData4Chris/Medium/abn10784_20111130_082334_3409.mat"
            #  != "../testdata/IIICData4Chris/Medium/sid1489_20140905_074518_10912.mat"
        ):
            continue

        # TODO: Deal with file loading issues
        data = iiic_detector.utils.MATFile(fname, eng)

        P_model = data.P_model
        Y_human = data.Y_human

        Y_human_numvotes = Y_human.iloc[0, :].sum()
        Y_human_p = Y_human / Y_human_numvotes

        Y_human_s = Y_human_p[Y_human_p >= cutoff_pct].dropna(axis=1)

        theseTypes = Y_human_s.columns.tolist()
        theseTypes_p = Y_human_s.iloc[0, :].tolist()

        thisResult = iiic_detector.utils.FileData(
            n=idx,
            file=fname,
            y_human_types=theseTypes,
            y_human_types_pct=theseTypes_p,
            stats_obj=None,
        )

        Y_human_s_max = Y_human_s.max(axis="columns").item()
        thisType = Y_human_s.idxmax(axis="columns").item()

        data_obj = None

        if not np.isnan(Y_human_s_max):
            if thisType in ["GRDA", "LRDA"]:
                result_table.append(thisResult)
            elif thisType in ["LPD", "GPD"]:
                # Remove EKG
                seg = np.asarray(data.data)[:19, :]
                freq, channels, interpretation, data_obj = iiic_detector.pd_detect(seg)
            elif thisType in ["Seizure"]:
                # Remove EKG
                seg = np.asarray(data.data)[:19, :]
                freq, channels, interpretation, data_obj = iiic_detector.pd_detect(seg)
            elif thisType in ["Other"]:
                result_table.append(thisResult)
            else:
                result_table.append(thisResult)

        # TODO: Generate Figure Section
        if isinstance(data_obj, dict):
            data_obj["num"] = idx
            data_obj["filename"] = fname
            data_obj["theseTypes"] = theseTypes
            data_obj["theseTypes_prctVotes"] = theseTypes_p
            data_obj["Y_human"] = Y_human
            data_obj["P_model"] = P_model

            f1, this_stats_obj = iiic_detector.show_eeg_events_and_stats(data_obj)
        pass
    return


if __name__ == "__main__":
    fire.Fire(main)
