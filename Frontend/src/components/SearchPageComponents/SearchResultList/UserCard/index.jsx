import React from "react";

import SearchCard from "../SearchCard";

// Card to display search result for a single user
const UserCard = ({ item: user }) => {
  return (
    <SearchCard
      avatar={user.image}
      title={user.name}
      secondary={` ${user.affiliation}`}
      tertiary={`${user.email}`}
    />
  );
};

export default UserCard;
