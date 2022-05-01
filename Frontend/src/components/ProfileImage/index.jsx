import React, { useState, useEffect } from 'react';
import Avatar from '@mui/material/Avatar';

import ProfileImageUploadDialog from 'components/general/upload/ProfileImageUploadDialog';
import styles from './profileImage.module.css';

import stringToColor from 'shared/utils/stringColor';
import getProfile from 'shared/services/profile';
import { useAuth } from 'shared/context/authContext';

function ProfileImage({ sizePixels, editable = false, overrideProfileImage = null }) {
  const [user, setUser] = useAuth();
  async function loadUser() {
    setUser(await getProfile());
  }
  useEffect(loadUser, [setUser]);

  const [uploadDialogOpen, setUploadDialogOpen] = useState(false);

  return (
    <>
      <div
        className={styles.profileImage} // TODO: load actual image
        style={{
          height: `${sizePixels}px`,
          width: `${sizePixels}px`,
        }}
      >
        <Avatar
          src={overrideProfileImage || user.avatarUrl}
          alt={`${user.firstName} ${user.lastName}`}
          sx={{
            backgroundColor: stringToColor(`${user.firstName} ${user.lastName}`),
            width: sizePixels,
            height: sizePixels,
          }}
        >
          {(user.firstName || '?')[0]}
        </Avatar>
        {editable && (
          <div
            className={styles.editButton}
            onClick={() => setUploadDialogOpen(true)}
            role="button"
            onKeyPress={() => setUploadDialogOpen(true)}
            tabIndex={0}
          >
            <span>Edit</span>
          </div>
        )}
      </div>
      {uploadDialogOpen && (
        <ProfileImageUploadDialog
          open={uploadDialogOpen}
          onClose={() => {
            setUploadDialogOpen(false);
            loadUser();
          }}
          onChange={(newUrl) => {
            // update without reload
            user.avatarUrl = newUrl;
            setUser(user);
          }}
        />
      )}
    </>
  );
}

export default ProfileImage;
