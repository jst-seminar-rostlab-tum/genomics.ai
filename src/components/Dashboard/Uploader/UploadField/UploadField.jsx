import React from 'react';
import styles from './uploadfield.module.css';
import UploadButton from '../UploadButton/UploadButton';

/*
The disabled prop makes the upload field disabled
*/

function UploadField({ onChange, disabled, selectedFile }) {
  return (
    <div className={styles.flexContainer}>
      <div className={styles.uploadButtonContainer}>
        <UploadButton
          onChange={onChange}
          selectedFile={selectedFile}
          disabled={disabled}
        />
      </div>
    </div>
  );
}

export default UploadField;
