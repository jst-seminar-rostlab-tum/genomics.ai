import React from 'react';

import SearchCard from '../SearchCard';

// Card to display search result for a single user
function UserCard({ item: user }) {
  return (
    <SearchCard
      avatar={user.avatarUrl}
      title={`${user.firstName} ${user.lastName}`}
      link={`/sequencer/users/${user._id}`}
      secondary={` ${user.note}`}
      displayAvatar
      // tertiary={`${user.email}`}
    />
  );
}

export default UserCard;
