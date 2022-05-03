import React, { useState } from 'react';
import Avatar from '@mui/material/Avatar';

import InstitutionProfileImageUploadDialog from 'components/general/upload/InstitutionProfileImageUploadDialog';
import styles from './institutionAvatar.module.css';

import stringToColor from 'shared/utils/stringColor';

function InstitutionAvatar({ editable = false, institution, onChange }) {
  const sizePixels = 200;

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
          src={institution.avatarUrl}
          alt={`${institution.name}`}
          sx={{
            backgroundColor: stringToColor(institution.name),
            width: sizePixels,
            height: sizePixels,
          }}
        >
          {(institution.name || '?')[0]}
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
        <InstitutionProfileImageUploadDialog
          institution={institution}
          open={uploadDialogOpen}
          onClose={() => {
            setUploadDialogOpen(false);
          }}
          onChange={onChange}
        />
      )}
    </>
  );
}

export default InstitutionAvatar;
