import React from "react";

import SearchCard from "../SearchCard";

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
