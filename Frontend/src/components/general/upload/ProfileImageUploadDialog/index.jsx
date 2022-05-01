import React, { useState, useEffect } from 'react';
import MemberCard from 'components/members/MemberCard';
import ImageUploadDialog from '../ImageUploadDialog';

import getProfile from 'shared/services/profile';

import { BACKEND_ADDRESS } from 'shared/utils/common/constants';
import { getAuthAndJsonHeader } from 'shared/utils/common/utils';

export default function ProfileImageUploadDialog({
  open, onClose, onChange,
}) {
  const [user, setUser] = useState({});
  useEffect(() => {
    getProfile()
      .then(setUser);
  });

  async function upload(blob) {
    const response = await fetch(`${BACKEND_ADDRESS}/user-avatar`, {
      method: 'POST',
      headers: {
        ...getAuthAndJsonHeader(),
        'content-type': 'image/png',
      },
      body: blob,
    });
    if (response.status !== 200) {
      if (response.status === 413) {
        // should not happen because of the max file size in the filedrop component,
        // but who knows, maybe the backend will change the limit in the future
        alert('This image is too large. Please provide one with a smaller file size.');
        return;
      }
      throw Error("Couldn't upload user avatar");
    }
    // response is the new profile image URL
    if (onChange) {
      onChange(await response.text());
    }
    onClose();
  }

  async function reset() {
    await fetch(`${BACKEND_ADDRESS}/user-avatar`, {
      method: 'DELETE',
      headers: getAuthAndJsonHeader(),
    });
    onClose();
  }

  return (
    <ImageUploadDialog
      open={open}
      onClose={onClose}
      title="Upload Profile Image"
      maxFileSizeMB={1}
      croppable
      preview={(imgURL) => (
        <MemberCard memberId={user.id} overrideProfilePicture={imgURL} />
      )}
      additionalButtons={[
        { text: 'Reset', func: reset },
      ]}
      onUpload={(blob) => upload(blob)}
    />
  );
}
