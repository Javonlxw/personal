/* Start Header ************************************************************************/
/*!
\file       Input_Handler.cpp
\author     William Yoong
\par email: william.yoong\@digipen.edu
\brief
Input from user handled here.

Copyright (C) 20xx DigiPen Institute of Technology.
Reproduction or disclosure of this file or its contents
without the prior written consent of DigiPen Institute of
Technology is prohibited.
*/
/* End Header **************************************************************************/

#include "Input_Handler.h"
#include "Camera.h"
#include <iostream>

/*

struct Input_Queue
{
	int value;
	Input_Queue *next;
};

// input buffer
namespace
{
	Input_Queue *user_input = nullptr;
	Input_Queue *input_tail = nullptr;
}

// create a new input
Input_Queue *New_Input(int value)
{
	Input_Queue *new_node = new Input_Queue;
	new_node->value = value;
	new_node->next = nullptr;

	return new_node;
}

// checks if the input buffer is empty or not
bool IsEmpty(void)
{
	return (user_input == nullptr);
}


// adds a new node to the end of the input buffer
void Enqueue(int value)
{
	// input buffer is empty
	if (IsEmpty())
	{
		user_input = New_Input(value);
		input_tail = user_input;
		return;
	}
	else // adds the node to the end of the input buffer
	{
		input_tail->next = New_Input(value);
		input_tail = input_tail->next;
	}


}


// removes the front of the input buffer
void Dequeue(void)
{
	// input buffer is empty
	if (IsEmpty())
		return;

	if (user_input->next == nullptr) // last item in input buffer
	{
		delete user_input;
		user_input = nullptr;
		return;
	}

	else 
	{
		Input_Queue *temp = user_input;
		user_input = temp->next;
		delete temp;
		return;
	}
}

// clears the input buffer/queue
void Clear_Buffer(void)
{
	if (IsEmpty())
		return;

	Input_Queue *temp = user_input;

	while (user_input)
	{
		user_input = user_input->next;
		delete temp;
		temp = user_input;
	}
	input_tail = nullptr;
}
*/
namespace Input {

	std::vector <int> input;
	
	const std::vector <int> &Get_User_Input(void)
	{
		using namespace Input;

		AEInputUpdate();

		if (AEInputCheckCurr(go_right))
		{
			input.push_back(go_right);
		}

		if (AEInputCheckCurr(go_left))
		{
			input.push_back(go_left);
		}

		if (AEInputCheckCurr(go_up))
		{
			input.push_back(go_up);
		}

		if (AEInputCheckCurr(go_down))
		{
			input.push_back(go_down);
		}

		if (AEInputCheckTriggered(jump_up))
		{
			input.push_back(jump_up);
		}

		if (AEInputCheckTriggered(fire))
		{
			input.push_back(fire);
		}

		if (AEInputCheckTriggered(special))
		{
			input.push_back(special);
		}

		// if esc key or the close button was triggered, stop the game
		if (AEInputCheckTriggered(quit_game) || 0 == AESysDoesWindowExist())
		{
			GSM::current = GS_QUIT;
		}

		return input;
	}

	void ClearBuffer(void)
	{
		input.clear();
	}
}