import React from 'react';
import profileDefault from 'assets/user.png';

import styles from './profileImage.module.css';

function ProfileImage({ sizePixels, editable = false }) {
  return (
    <div
      className={styles.profileImage} // TODO: load actual image
      style={{
        height: `${sizePixels}px`,
        width: `${sizePixels}px`,
        backgroundImage: `url(${profileDefault})`,
      }}
    >
      {editable && (
        <div className={styles.editButton}>
          <span>Edit</span>
        </div>
      )}
    </div>
  );
}

export default ProfileImage;
