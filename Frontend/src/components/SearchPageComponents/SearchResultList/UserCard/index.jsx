import React from "react";

import SearchCard from "../SearchCard";

// Card to display search result for a single user
const UserCard = ({ item: user }) => {
  return (
    <SearchCard
      avatar={user.image}
      title={user.name}
      link={`/sequencer/users/${user.id}`}
      secondary={` ${user.affiliation}`}
      tertiary={`${user.email}`}
    />
  );
};

export default UserCard;
