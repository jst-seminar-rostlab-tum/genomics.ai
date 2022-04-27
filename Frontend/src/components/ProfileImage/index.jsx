import React, { useState } from 'react';
import profileDefault from 'assets/user.png';

import ProfileImageUploadDialog from 'components/general/upload/ProfileImageUploadDialog';

import styles from './profileImage.module.css';

function ProfileImage({ sizePixels, editable = false, overrideProfileImage = null }) {
  const [uploadDialogOpen, setUploadDialogOpen] = useState(false);

  const openUploadDialog = () => setUploadDialogOpen(true);

  return (
    <>
      <div
        className={styles.profileImage} // TODO: load actual image
        style={{
          height: `${sizePixels}px`,
          width: `${sizePixels}px`,
          backgroundImage: `url(${overrideProfileImage || profileDefault})`,
        }}
      >
        {editable && (
          <div
            className={styles.editButton}
            onClick={() => openUploadDialog()}
            role="button"
            onKeyPress={() => openUploadDialog()}
            tabIndex={0}
          >
            <span>Edit</span>
          </div>
        )}
      </div>
      {uploadDialogOpen && (
        <ProfileImageUploadDialog
          open={uploadDialogOpen}
          onClose={() => setUploadDialogOpen(false)}
        />
      )}
    </>
  );
}

export default ProfileImage;
