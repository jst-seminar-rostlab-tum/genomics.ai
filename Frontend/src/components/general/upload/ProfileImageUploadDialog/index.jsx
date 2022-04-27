import React, { useState, useEffect } from 'react';
import MemberCard from 'components/members/MemberCard';
import ImageUploadDialog from '../ImageUploadDialog';

import getUser from 'shared/services/mock/user';

export default function ProfileImageUploadDialog({ open, onClose }) {
  const [user, setUser] = useState({});
  useEffect(() => {
    getUser()
      .then(setUser);
  });

  return (
    <ImageUploadDialog
      open={open}
      onClose={onClose}
      title="Upload Profile Image"
      maxFileSizeMB={5}
      croppable
      preview={(imgURL) => (
        <MemberCard memberId={user.id} overrideProfilePicture={imgURL} />
      )}
    />
  );
}
