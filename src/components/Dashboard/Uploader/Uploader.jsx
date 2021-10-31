import {
  Tooltip,
} from '@mui/material';
import React from 'react';

import HelpIcon from '@mui/icons-material/Help';
import styles from './uploader.module.css';
import SubmitForm from './SubmitForm/SubmitForm';

function Uploader() {
  return (
    <div className={styles.Container}>
      <div className={styles.topBar}>
        <div className={styles.title}>
          <h1>Upload your files</h1>
        </div>
        <div className={styles.tooltip}>
          <Tooltip title="Accepted File Format: .txt" placement="bottom-end">
            <HelpIcon />
          </Tooltip>
        </div>
      </div>
      <SubmitForm />
    </div>
  );
}

export default Uploader;
