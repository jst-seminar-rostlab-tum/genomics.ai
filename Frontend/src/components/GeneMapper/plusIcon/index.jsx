import { IconButton } from '@mui/material';
import Icon from '@mui/material/Icon';
import React from 'react';
import { useHistory } from 'react-router-dom';

export default function PlusIcon() {
  const history = useHistory();
  //  const handleClickWeb = () => {
  //    window.open('http://google.com', '_blank');

  return (
    // change the size of it
    // route to the right page!!
    <IconButton aria-label="plus icon" onClick={() => { history.push('./genemapper/selection'); }}>
      <Icon className="plusicon" color="primary" sx={{ fontSize: 25 }}>add_circle</Icon>
    </IconButton>

  );
}
